#pragma once

#include <eigen3/Eigen/Core>

#include <vector>
#include <memory>
#include <unordered_set>
#include <queue>
#include <algorithm>

#include "planning_algorithm/GridNode.h"

namespace planning_algorithm
{

auto greater_gridnode = [] (GridNodePtr a, GridNodePtr b) { return (a->fscore_ > b->fscore_); };

static const std::vector<Eigen::Vector3i> CONNECTION = 
{
    {-1, -1, 1},
    {-1, 0, 1},
    {-1, 1, 1},
    {0, 1, 1},
    {1, 1, 1},
    {1, 0, 1},
    {1, -1, 1},
    {0, -1, 1},
    {0, 0, 1},
    {-1, -1, 0},
    {-1, 0, 0},
    {-1, 1, 0},
    {0, 1, 0},
    {1, 1, 0},
    {1, 0, 0},
    {1, -1, 0},
    {0, -1, 0},
    {-1, -1, -1},
    {-1, 0, -1},
    {-1, 1, -1},
    {0, 1, -1},
    {1, 1, -1},
    {1, 0, -1},
    {1, -1, -1},
    {0, -1, -1},
    {0, 0, -1}
};

class Astar
{
public:
    
    Astar(double resolution, const Eigen::Vector3d& global_xyz_low, const Eigen::Vector3d& global_xyz_upp);
    ~Astar();

    void SetObstacle(const Eigen::Vector3i& idx);
    void SetObstacle(const Eigen::Vector3d& pt);
    bool SetStart(const Eigen::Vector3i& index);
    bool SetStart(const Eigen::Vector3d& coord);
    bool SetGoal(const Eigen::Vector3i& index);
    bool SetGoal(const Eigen::Vector3d& coord);
    bool GetPath(std::vector<Eigen::Vector3d>& path);

private:
    GridNodePtr*** grid_map_;
    GridNodePtr start_;
    GridNodePtr goal_;
    int GRID_XSIZE, GRID_YSIZE, GRID_ZSIZE;
    int GRID_YZSIZE, GRID_XYZSIZE;
    double GLOBAL_XLOW, GLOBAL_YLOW, GLOBAL_ZLOW;
    double GLOBAL_XUPP, GLOBAL_YUPP, GLOBAL_ZUPP;
    double resolution_, inv_resolution_;
    std::vector<Eigen::Vector3i> connection_;

    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, decltype(greater_gridnode)> openset_;
    std::unordered_set<GridNodePtr> closedset_;

    bool ReconstructPath(std::vector<Eigen::Vector3d>& path, GridNodePtr current);
    void GetNeighborNodes(std::vector<GridNodePtr>& neighbor_nodes, GridNodePtr current);
    double GetHeuristic(GridNodePtr current);
    Eigen::Vector3d Index2Coord(const Eigen::Vector3i& idx);
    Eigen::Vector3i Coord2Index(const Eigen::Vector3d& pt);
    
    bool Occupied(int x, int y, int z);
    bool Occupied(GridNodePtr current);
};

} /* namespace */