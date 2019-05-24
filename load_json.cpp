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

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>

#include "common.hpp"

degree_distribution load_degree_distribution(std::string file_name){
    boost::property_tree::ptree pt;
    uint64_t x;
    double y;
    std::vector<distribution_pair> lambda;
    std::vector<distribution_pair> rho;
    distribution_pair dpair;
    degree_distribution rtnv;

    read_json(file_name, pt);

    // extract variable node degree distribution
    BOOST_FOREACH(const boost::property_tree::ptree::value_type& child, pt.get_child("degree_distribution.lambda")) {
        const boost::property_tree::ptree& info = child.second;

        if(boost::optional<uint64_t> degree = info.get_optional<uint64_t>("degree")) {
            x = degree.get();
            std::cout << "degree: " << x << std::endl;
        } else {
            std::cerr << "degree is not found." << std::endl;
            exit(-1);
        }

        if(boost::optional<double> weight = info.get_optional<double>("weight")) {
            y = weight.get();
            std::cout << "weight: " << y << std::endl;
        } else {
            std::cerr << "weight is not found." << std::endl;
            exit(-1);
        }
        dpair = std::make_pair(x, y);
        lambda.push_back(dpair);
    }

    // extract check node degree distribution
    BOOST_FOREACH(const boost::property_tree::ptree::value_type& child, pt.get_child("degree_distribution.rho")) {
        const boost::property_tree::ptree& info = child.second;

        if(boost::optional<uint64_t> degree = info.get_optional<uint64_t>("degree")) {
            x = degree.get();
            std::cout << "degree: " << x << std::endl;
        } else {
            std::cerr << "degree is not found." << std::endl;
            exit(-1);
        }

        if(boost::optional<double> weight = info.get_optional<double>("weight")) {
            y = weight.get();
            std::cout << "weight: " << y << std::endl;
        } else {
            std::cerr << "weight is not found." << std::endl;
            exit(-1);
        }
        dpair = std::make_pair(x, y);
        rho.push_back(dpair);
    }

    rtnv = std::make_pair(lambda, rho);
    return(rtnv);
}

int main(int argc, char* argv[]) {
    std::vector<distribution_pair> pairv;
    std::vector<distribution_pair>::iterator vit;
    degree_distribution dist;

    if(argc != 2) {
        std::cerr << "Usage ./load_json <json file name>" << std::endl;
        return(EXIT_FAILURE);
    }
    dist = load_degree_distribution(std::string(argv[1]));

    // print lambda
    pairv = dist.first;
    for(vit = pairv.begin(); vit != pairv.end(); ++vit) {
        std::cout << "degree = " << vit->first << "\tweght = " << vit->second << std::endl;
    }
    // print rho
    pairv = dist.second;
    for(vit = pairv.begin(); vit != pairv.end(); ++vit) {
        std::cout << "degree = " << vit->first << "\tweght = " << vit->second << std::endl;
    }
    return(EXIT_SUCCESS);
}