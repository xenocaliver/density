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
#ifndef PLOT_H
#define PLOT_H

#include <string>
#include <QwtPlot>

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