#ifndef PLOT_H
#define PLOT_H

#include <string>
#include <qwt_plot.h>

class QwtPlotCurve;

class Plot : public QwtPlot
{
    Q_OBJECT
public:
    explicit Plot(QWidget *parent = 0);
    void plotCurve(double*, double*, uint64_t);
    void setupCurve(std::string);
    void updateCurve(double*, double*, uint64_t);

signals:
    void emitSignal(void);

private:
    QwtPlotCurve *curve1_;
};

#endif // PLOT_H