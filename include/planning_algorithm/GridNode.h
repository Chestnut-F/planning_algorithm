#pragma once

#include <eigen3/Eigen/Core>
#include <memory>

namespace planning_algorithm
{

static const double INF = 1e20;

class GridNode;
typedef std::shared_ptr<GridNode> GridNodePtr;

class GridNode
{
public:
    GridNode(const Eigen::Vector3i& index, const Eigen::Vector3d& coord):
        index_(index), coord_(coord),
        gscore_(INF), hscore_(INF), fscore_(INF),
        traversed_(false), occupied_(false)
    {
        comefrom_ = nullptr;
    }

    Eigen::Vector3i index_;
    Eigen::Vector3d coord_;

    GridNodePtr comefrom_;

    double gscore_, hscore_, fscore_;

    bool traversed_;
    bool occupied_;
};



} /* namespace */