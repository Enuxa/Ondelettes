#ifndef WAVELETS_H
#define WAVELETS_H

#include <complex>
#include <QVector>
#include <functional>
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
    int k;
    bool operator<(const Coord &coord) const {
        if(n < coord.n)
            return true;
        if(k < coord.k)
            return true;
        return false;
    }
};

typedef double Vector;
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
public:

    virtual Complex psi(double t) = 0;
    virtual Complex psiPrimitive(double x) = 0;
    virtual Complex phi(double t) = 0;
    virtual Complex phiPrimitive(double x) = 0;
    virtual Complex1DFunc analyzingFunction(int n, int k) {
        Complex1DFunc f = [=] (double t) {
            return pow(2, -n/2.0)*psi(pow(2,-n)*t-k);
        };
        return f;
    }

    virtual Coefficients coefficientsOf(Complex2DFunc f) {
        Coefficients coeffs;
        for(int n=1;n<=9;n++) {
            for(int k=0;k<pow(2, n)*N;k++) {
                Coord coord;
                coord.n = n;
                coord.k = k;
                Complex coeff = std::polar(0.0);
                for(int x=0; x<N; x++) {
                    double t11 = pow(2, -n)*x-k;
                    double t22 = pow(2, -n)*(x+1)-k;
                    Complex z1 = psiPrimitive(t11);
                    Complex z2 = psiPrimitive(t22);
                    Complex fx = f(x);
                    Complex coeffp = pow(2, 0.5*n)*fx*(z2-z1);;
                    coeff += coeffp;
                }
                coeffs[coord] = coeff;
            }
        }
        Coord coord;
        coord.n = -1;
        coord.k = 0;
        Complex coeff = std::polar(0.0);
        int n=9,k=0;
        for(int x=0; x<N; x++) {
            double t11 = pow(2, -n)*x-k;
            double t22 = pow(2, -n)*(x+1)-k;
            Complex z1 = phiPrimitive(t11);
            Complex z2 = phiPrimitive(t22);
            Complex fx = f(x);
            Complex coeffp = pow(2, 0.5*n)*fx*(z2-z1);;
            coeff += coeffp;
        }
        coeffs[coord] = coeff;
        return coeffs;
    }

    virtual Complex2DFunc reconstructFunction(Coefficients coeffs) {
        Complex2DFunc f = [=] (Vector x) {
            Complex value = std::polar(0.0);
            foreach(Coord coord, coeffs.keys()) {
                if(coord.n == -1)
                    value += coeffs[coord]*pow(2, -0.5*9)*phi(pow(2, -9)*x);
                else
                    value += coeffs[coord]*analyzingFunction(coord.n,coord.k)(x);
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
    virtual Complex phi(double x) {
        if(x>=0 && x<1)
            return std::polar(1.0);
        return std::polar(0.0);
    }
    virtual Complex phiPrimitive(double x) {
        if(x>=0 && x<1)
            return std::polar(x);
        else if(x>=1)
            return std::polar(1.0);
        return std::polar(0.0);
    }
};

template <int N>
class MexicanHatBase : public WaveletsBase<N> {
public:
    double lambda = pow(2,0.75)/(sqrt(3)*pow(PI,0.25));
    virtual Complex psi(double t) {
        //if(t<0 || t>=1) return std::polar(0.0);
        return std::polar(lambda*exp(-t*t*0.5)*(1-t*t));
    }
    virtual Complex psiPrimitive(double t) {
        //if(t<0 || t>=1) return std::polar(0.0);
        return std::polar(lambda*exp(-t*t*0.5)*t);
    }
    virtual Complex phi(double x) {
        if(x>=0 && x<1)
            return std::polar(1.0);
        return std::polar(0.0);
    }
    virtual Complex phiPrimitive(double x) {
        if(x>=0 && x<1)
            return std::polar(x);
        else if(x>=1)
            return std::polar(1.0);
        return std::polar(0.0);
    }
};


#endif // WAVELETS_H
