#pragma once

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <vector>
#include <utility>

namespace lib_rtree
{
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;

    // Boost.Geometry用の型定義
    typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
    typedef bg::model::box<point_t> box_t;
    
    // R-treeに格納する値: bounding boxと線分の端点のペア
    typedef std::pair<box_t, std::pair<std::vector<double>, std::vector<double>>> rtree_value_t;
}

