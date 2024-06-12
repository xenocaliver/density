/* 
* This file is part of the density distribution (https://github.com/xenocaliver/density).
* Copyright (c) 2019 Akiyoshi Hashimoto.
* 
* This program is free software: you can redistribute it and/or modify  
* it under the terms of the GNU General Public License as published by  
* the Free Software Foundation, version 3.
*
* This program is distributed in the hope that it will be useful, but 
* WITHOUT ANY WARRANTY; without even the implied warranty of 
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License 
* along with this program. If not, see <http://www.gnu.org/licenses/>
*/
#include "plot.hpp"
#include <qwt_plot_curve.h>
#include <qwt_legend.h>
#include <qwt_plot_grid.h>

Plot::Plot(QWidget *parent) :
    QwtPlot(parent)
{
    setCanvasBackground(QColor(Qt::white));
    // legend
    insertLegend(new QwtLegend(), QwtPlot::BottomLegend);

    // grid 
    QwtPlotGrid *grid1 = new QwtPlotGrid;
    grid1->enableXMin(true);
    grid1->setMajorPen(Qt::black, 0, Qt::DashDotLine);
    grid1->setMinorPen(Qt::gray, 0 , Qt::DotLine);
    grid1->attach(this);

    // axis 
    setAxisTitle(QwtPlot::xBottom, "x");
    setAxisTitle(QwtPlot::yLeft, "y");

}

void Plot::plotCurve(double* x, double* y, uint64_t array_size)
{
    curve1_->setSamples(x, y, array_size);
}

void Plot::setupCurve(std::string CurveName) {
    // setup curve
    QString str = QString::fromStdString(CurveName);
    curve1_ = new QwtPlotCurve(str);
    curve1_->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve1_->setPen(Qt::red);
    curve1_->attach(this);
}

void Plot::updateCurve(double* x, double* y, uint64_t array_size) {
    curve1_->setSamples(x, y, array_size);
}
#include "moc_plot.cpp"