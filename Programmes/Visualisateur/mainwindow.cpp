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

    ComplexFunc ftest = [=] (double t) {
        return std::polar(cos(2*PI*12*t));
    };

    FourierBase base;
    Coordinates2D coords = base.coordinatesOf(ftest, 50, 2);

    /*ui->plotUp->insertLegend( new QwtLegend() );
    QwtPlotGrid *grid = new QwtPlotGrid();
    grid->attach(ui->plotUp);*/

    QwtPlotCurve *curve = new QwtPlotCurve();
    curve->setPen(Qt::blue, 4);
    curve->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    QwtSymbol *symbol = new QwtSymbol( QwtSymbol::Ellipse,
        QBrush( Qt::yellow ), QPen( Qt::red, 2 ), QSize( 8, 8 ) );
    curve->setSymbol(symbol);

    QPolygonF points;
    for(int i=0;i<N;i++) {
        double t = ((double)i)/N;
        points << QPointF(t, ftest(t).real());
    }
    curve->setSamples(points);

    curve->attach(ui->plotUp);

    QwtMatrixRasterData *rasterData = new QwtMatrixRasterData();
    for(int a=0;a<50;a++) {
        for(int b=0;b<2;b++) {
            rasterData->setValue(b, a, coords[a][b].real());
        }
    }

    QwtLinearColorMap *colorMap = new QwtLinearColorMap(Qt::darkCyan, Qt::red);
    colorMap->addColorStop(0.0, Qt::cyan);
    colorMap->addColorStop(0.6, Qt::green);
    colorMap->addColorStop(1.0, Qt::yellow);

    QwtPlotSpectrogram *spectrogram = new QwtPlotSpectrogram();
    spectrogram->setColorMap(colorMap);
    spectrogram->setData(rasterData);
    //spectrogram->setDisplayMode(QwtPlotSpectrogram::ImageMode, true);

    spectrogram->attach(ui->plotDown);

}

MainWindow::~MainWindow()
{
    delete ui;
}
