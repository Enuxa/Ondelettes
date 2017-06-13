#ifndef WAVELETS_H
#define WAVELETS_H

#include <complex>
#include <QVector>
#include <functional>
#include <QGenericMatrix>
#include <QList>
#include <QMap>
#include <QtDebug>

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

class Coord {
public:
    int n;
    QGenericMatrix<1,2,double> k;
    int i;
    bool operator<(const Coord &coord) const {
        if(n < coord.n)
            return true;
        if(k(0,0) < coord.k(0,0))
            return true;
        if(k(1,0) < coord.k(1,0))
            return true;
        if(i < coord.i)
            return true;
        return false;
    }
};

typedef QGenericMatrix<1,2,double> Vector;
typedef std::function<Complex(double)> Complex1DFunc; // Type ComplexFunc représentant une fonction f:double->Complex
typedef std::function<Complex(Vector)> Complex2DFunc; // Type ComplexFunc représentant une fonction f:Vector->Complex
typedef QMap<Coord, Complex> Coefficients; // Type représentant une liste de Complex

/**
 * @brief Produit scalaire entre deux fonctions
 */
Complex scalarProduct(Complex1DFunc f, Complex1DFunc g, int N) {
    Complex value;
    for(int i=0;i<N;i++) {
        value += (f(i+0.0)*std::conj(g(i+0.0)) + f(i+1.0)*std::conj(g(i+1.0)));
    }
    return value;
}


/**
 * @brief Classe abstraite dont les classes filles representent une base de Hilbert de L_2([0,N]^2)
 */
template <int N>
class WaveletsBase {
protected:
    QList<QGenericMatrix<1,2,double>> translationVectors;
public:
    WaveletsBase() {
        translationVectors.append(Vector(new double[2]{1.0,0.0}));
        translationVectors.append(Vector(new double[2]{1.0,1.0}));
        translationVectors.append(Vector(new double[2]{0.0,1.0}));
    }

    virtual Complex psi(double t) = 0;
    //virtual Complex psiPrimitive(double x) = 0;
    virtual Complex1DFunc analyzingFunction(int n) {
        Complex1DFunc f = [=] (double t) {
            return pow(2, n)*psi(pow(2,n)*t);
        };
        return f;
    }

    virtual Coefficients coefficientsOf(Complex2DFunc f) {
        Coefficients coeffs;
        for(int n=0;n<5;n++) {
            for(int ka=0;ka<N;ka++) {
                for(int kb=0;kb<N;kb++) {
                    for(int i=0;i<3;i++) {
                        Coord coord;
                        coord.n = n;
                        coord.i = i;
                        coord.k = Vector(new double[2]{(double)ka,(double)kb});
                        Vector tv = translationVectors[i];
                        Complex1DFunc fRotated = [=] (double t) {
                            return f(coord.k+t*tv);
                        };
                        Complex coeff = scalarProduct(fRotated, analyzingFunction(n), N);
                        coeffs[coord] = coeff;

                        qDebug() << coord.n << (int)(coeff*256);
                    }
                }
            }
        }
        return coeffs;
    }

    virtual Complex2DFunc reconstructFunction(Coefficients coeffs) {
        Complex2DFunc f = [=] (Vector x) {
            Complex value = std::polar(0.0);
            foreach(Coord coord, coeffs.keys()) {
                //value += coeffs[coord]*analyzingFunction(coord.n)();
            }
            return value;
        };
        return f;
    }
};

template <int N>
class HaarBase : public WaveletsBase<N> {
public:
    virtual Complex psi(double x) {
        if(x>=0 && x<0.5)
            return std::polar(1.0);
        if(x>=0.5 && x<1)
            return std::polar(-1.0);
        return std::polar(0.0);
    }
    virtual Complex psiPrimitive(double x) {
        if(x>=0 && x<0.5)
            return std::polar(x);
        if(x>=0.5 && x<1)
            return std::polar(1.0-x);
        return std::polar(0.0);
    }
};


#endif // WAVELETS_H
