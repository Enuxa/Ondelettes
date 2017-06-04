#ifndef WAVELETS_H
#define WAVELETS_H

#include <complex>
#include <QVector>
#include <functional>
#include <QtDebug>

#define N 100
#ifndef PI
   #define PI 3.14159265358979323846
#endif

/**
 * @brief Classe représentant un nombre complexe
 */
class Complex : public std::complex<double> {
public:
    Complex()
        : std::complex<double>(0.0) {
    }
    Complex(std::complex<double> c)
        : std::complex<double>(c) {
    }
};

/**
 * @brief Fonction servant à afficher un nombre complexe
 */
QDebug operator<<(QDebug debug, const Complex c)
{
    QDebugStateSaver saver(debug);
    debug.nospace() <<  c.real() << "+" << c.imag() << 'i';

    return debug;
}

//typedef Complex(*ComplexFunc)(double);
typedef std::function<Complex(double)> ComplexFunc; // Type ComplexFunc représentant une fonction f:double->Complex
typedef QVector<Complex> Coordinates1D; // Type Coordinates1D représentant une liste de Complex
typedef QVector<Coordinates1D> Coordinates2D; // Type Coordinates1D représentant une liste de liste de Complex

/**
 * @brief Produit scalaire entre deux fonctions
 */
Complex scalarProduct(ComplexFunc f, ComplexFunc g) {
    Complex value;
    for(int i=0;i<N;i++) {
        double t = ((double)i)/N;
        value += (f(t)*std::conj(g(t)) + f(t+1)*std::conj(g(t+1)))/std::polar(2.0*N);
    }
    return value;
}

/**
 * @brief Classe abstraite dont les classes filles representent une base de Hilbert de L_2([0,1])
 */
class HilbertBase {
public:
    /**
     * @brief Cette fonction attend deux arguments a et b (position et dilatation),
     *        et retourne la fonction analysante correspondante
     */
    virtual ComplexFunc analyzingFunction(int a, int b) = 0;

    /**
     * @brief Retourne le tableau des coordonnees de f par rapport a la base tronquee de (0,0) a (A-1,B-1)
     */
    virtual Coordinates2D coordinatesOf(ComplexFunc f, int A, int B) {
        Coordinates2D coordsMatrix;
        for(int a=0;a<A;a++) {
            Coordinates1D coordsLine;
            for(int b=0;b<B;b++) {
                coordsLine.append(scalarProduct(f, analyzingFunction(a, b)));
            }
            coordsMatrix.append(coordsLine);
        }
        return coordsMatrix;
    }

    /**
     * @brief Reconstruit une fonction à partir de ses coordonnées
     */
    virtual ComplexFunc reconstructFunction(Coordinates2D coordsMatrix) {
        ComplexFunc f = [=] (double t) {
            Complex value = std::polar(0.0);
            int A = coordsMatrix.size();
            for(int a=0;a<A;a++) {
                int B = coordsMatrix[a].size();
                for(int b=0;b<B;b++) {
                    value += coordsMatrix[a][b]*analyzingFunction(a, b)(t);
                }
            }
            return value;
        };
        return f;
    }
};

class FourierBase : public HilbertBase {
public:
    virtual ComplexFunc analyzingFunction(int a, int b) {
        ComplexFunc f = [=] (double t) {
            if(b==0)
                return std::polar(1.0, 2*PI*a*t);
            else
                return std::polar(0.0, 0.0);
        };
        return f;
    }
};


#endif // WAVELETS_H
