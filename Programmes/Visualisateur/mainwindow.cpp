#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "wavelets.h"
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>
#include <qwt_symbol.h>
#include <qwt_legend.h>
#include <qwt_plot_rasteritem.h>
#include <qwt_matrix_raster_data.h>
#include <qwt_plot_spectrogram.h>
#include <qwt_color_map.h>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QImage imagetest(":/imagetest.pgm");

    qDebug() << imagetest.size();

    Complex1DFunc f = [=] (Vector x) {
        //return std::polar(cos(cos(x)*x));
        //return std::polar(cos(x*x/100)*exp(-x*x/1000));
        return std::polar((double)(imagetest.pixel(x,1000)%256));
    };

    MexicanHatBase<200> base;
    Coefficients coeffs = base.coefficientsOf(f);

    Complex1DFunc g = base.reconstructFunction(coeffs);

    QwtPlotCurve *curve = new QwtPlotCurve();
    curve->setPen(Qt::blue, 1);
    curve->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    /*QwtSymbol *symbol = new QwtSymbol( QwtSymbol::Ellipse,
        QBrush( Qt::yellow ), QPen( Qt::red, 2 ), QSize( 8, 8 ) );
    curve->setSymbol(symbol);*/

    QPolygonF points;
    for(int i=0;i<200;i++) {
        points << QPointF(i, f(i).real());
    }
    curve->setSamples(points);

    curve->attach(ui->plotUp);

    QwtPlotCurve *curve2 = new QwtPlotCurve();
    curve2->setPen(Qt::red, 1);
    curve2->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    /*QwtSymbol *symbol = new QwtSymbol( QwtSymbol::Ellipse,
        QBrush( Qt::yellow ), QPen( Qt::red, 2 ), QSize( 8, 8 ) );
    curve->setSymbol(symbol);*/

    QPolygonF points2;
    for(int i=0;i<200;i++) {
        points2 << QPointF(i, g(i).real());
    }
    curve2->setSamples(points2);

    curve2->attach(ui->plotUp);

}

MainWindow::~MainWindow()
{
    delete ui;
}
