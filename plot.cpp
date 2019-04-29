#include "plot.hpp"
#include <qwt_plot_curve.h>
#include <qwt_legend.h>
#include <qwt_plot_grid.h>

Plot::Plot(QWidget *parent) :
    QwtPlot(parent)
{
    setCanvasBackground(QColor(Qt::white));
    // 凡例
    insertLegend(new QwtLegend(), QwtPlot::BottomLegend);

    // グリッド
    QwtPlotGrid *grid1 = new QwtPlotGrid;
    grid1->enableXMin(true);
    grid1->setMajorPen(QPen(Qt::black, 0, Qt::DashDotLine));
    grid1->setMinorPen(QPen(Qt::gray, 0 , Qt::DotLine));
    grid1->attach(this);

    // 軸
    setAxisTitle(QwtPlot::xBottom, "x");
    setAxisTitle(QwtPlot::yLeft, "y");

}

void Plot::plotCurve(double* x, double* y, uint64_t array_size)
{
    curve1_->setSamples(x, y, array_size);
}

void Plot::setupCurve(std::string CurveName) {
    // 曲線の設定
    QString str = QString::fromStdString(CurveName);
    curve1_ = new QwtPlotCurve(str);
    curve1_->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve1_->setPen(QPen(Qt::red));
    curve1_->attach(this);
}

void Plot::updateCurve(double* x, double* y, uint64_t array_size) {
    curve1_->setSamples(x, y, array_size);
}
#include "moc_plot.cpp"