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
#ifndef QVIEWPDF_HPP
#define QVIEWPDF_HPP
#include <QtCharts/QChartView>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QLegendMarker>
#include <QtGui/QImage>
#include <QtGui/QPainter>

QT_CHARTS_USE_NAMESPACE

class Viewpdf : public QChartView {
    Q_OBJECT
public:
    Viewpdf(QWidget *parent = 0) : QChartView(new QChart(), parent) {
        QScatterSeries *series = new QScatterSeries();
    }
    ~Viewpdf(){}
};
#endif // QVIEWPDF_HPP