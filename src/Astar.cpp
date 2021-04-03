#include "planning_algorithm/Astar.h"

namespace planning_algorithm
{

Astar::Astar(double resolution, const Eigen::Vector3d& global_xyz_low, const Eigen::Vector3d& global_xyz_upp):
    openset_(greater_gridnode)
{
    GLOBAL_XLOW = global_xyz_low(0);
    GLOBAL_YLOW = global_xyz_low(1);
    GLOBAL_ZLOW = global_xyz_low(2);

    GLOBAL_XUPP = global_xyz_upp(0);
    GLOBAL_YUPP = global_xyz_upp(1);
    GLOBAL_ZUPP = global_xyz_upp(2);

    resolution_ = resolution;
    inv_resolution_ = 1.0 / resolution_;

    GRID_XSIZE = static_cast<int>( (GLOBAL_XUPP - GLOBAL_XLOW) * inv_resolution_ );
    GRID_YSIZE = static_cast<int>( (GLOBAL_YUPP - GLOBAL_YLOW) * inv_resolution_ );
    GRID_ZSIZE = static_cast<int>( (GLOBAL_ZUPP - GLOBAL_ZLOW) * inv_resolution_ );
    GRID_YZSIZE = GRID_YSIZE * GRID_ZSIZE;
    GRID_XYZSIZE = GRID_XSIZE * GRID_YZSIZE;

    grid_map_ = new GridNodePtr** [GRID_XSIZE];
    for (int i = 0; i < GRID_XSIZE; ++i)
    {
        grid_map_[i] = new GridNodePtr* [GRID_YSIZE];
        for (int j = 0; j < GRID_YSIZE; ++j)
        {
            grid_map_[i][j] = new GridNodePtr [GRID_ZSIZE];
            for (int k = 0; k < GRID_ZSIZE; ++k)
            {
                Eigen::Vector3i tmpIdx(i, j, k);
                Eigen::Vector3d pos = Index2Coord(tmpIdx);
                grid_map_[i][j][k] = std::make_shared<GridNode>(tmpIdx, pos);
            }
        }
    }
}

Astar::~Astar()
{

}

Eigen::Vector3d Astar::Index2Coord(const Eigen::Vector3i& idx)
{
    Eigen::Vector3d pt;

    pt(0) = ((double)idx(0) + 0.5) * resolution_ + GLOBAL_XLOW;
    pt(1) = ((double)idx(1) + 0.5) * resolution_ + GLOBAL_YLOW;
    pt(2) = ((double)idx(2) + 0.5) * resolution_ + GLOBAL_ZLOW;

    return pt;
}

Eigen::Vector3i Astar::Coord2Index(const Eigen::Vector3d& pt)
{
    Eigen::Vector3i idx;
    
    idx <<  std::min( std::max( int( (pt(0) - GLOBAL_XLOW) * inv_resolution_), 0), GRID_XSIZE - 1),
            std::min( std::max( int( (pt(1) - GLOBAL_YLOW) * inv_resolution_), 0), GRID_YSIZE - 1),
            std::min( std::max( int( (pt(2) - GLOBAL_ZLOW) * inv_resolution_), 0), GRID_ZSIZE - 1);

    return idx;
}

void Astar::SetObstacle(const Eigen::Vector3i& index)
{
    if (index(0) >= 0 && index(0) < GRID_XSIZE && 
        index(1) >= 0 && index(1) < GRID_YSIZE && 
        index(2) >= 0 && index(2) < GRID_ZSIZE)
    {
        grid_map_[index(0)][index(1)][index(2)]->occupied_ = true;
    }
    else
    {
        return;
    }
}

void Astar::SetObstacle(const Eigen::Vector3d& pt)
{
    if (pt(0) < GLOBAL_XUPP && pt(1) < GLOBAL_YUPP && pt(2) < GLOBAL_ZUPP && 
        pt(0) > GLOBAL_XLOW && pt(1) > GLOBAL_YLOW && pt(2) > GLOBAL_ZLOW)
    {
        Eigen::Vector3i idx = Coord2Index(pt);
        grid_map_[idx(0)][idx(1)][idx(2)]->occupied_ = true;
    }
    else
    {
        return;
    }
}

bool Astar::SetStart(const Eigen::Vector3i& index)
{
    if (index(0) >= 0 && index(0) < GRID_XSIZE && 
        index(1) >= 0 && index(1) < GRID_YSIZE && 
        index(2) >= 0 && index(2) < GRID_ZSIZE &&
        !Occupied(index(0), index(1), index(2)))
    {
        start_ = grid_map_[index(0)][index(1)][index(2)];
        return true;
    }
    else
    {
        return false;
    }
}

bool Astar::SetStart(const Eigen::Vector3d& coord)
{
    Eigen::Vector3i idx = Coord2Index(coord);
    return SetStart(idx);
}

bool Astar::SetGoal(const Eigen::Vector3i& index)
{
    if (index(0) >= 0 && index(0) < GRID_XSIZE && 
        index(1) >= 0 && index(1) < GRID_YSIZE && 
        index(2) >= 0 && index(2) < GRID_ZSIZE &&
        !Occupied(index(0), index(1), index(2)))
    {
        goal_ = grid_map_[index(0)][index(1)][index(2)];
        return true;
    }
    else
    {
        return false;
    }
}

bool Astar::SetGoal(const Eigen::Vector3d& coord)
{
    Eigen::Vector3i idx = Coord2Index(coord);
    return SetGoal(idx);
}

bool Astar::GetPath(std::vector<Eigen::Vector3d>& path)
{
    start_->gscore_ = 0.;
    start_->fscore_ = GetHeuristic(start_);
    openset_.push(start_);

    while (!openset_.empty())
    {
        GridNodePtr current = openset_.top();

        if (current == goal_)
            return ReconstructPath(path, current);
        
        openset_.pop();
        closedset_.insert(current);
        
        std::vector<GridNodePtr> neighbor_nodes;
        GetNeighborNodes(neighbor_nodes, current);

        for (GridNodePtr neighbor: neighbor_nodes)
        {
            if (closedset_.count(neighbor))
            {
                continue;
            }
            else
            {
                double tentative_gscore = current->gscore_ + (current->coord_ - neighbor->coord_).norm();
                if (tentative_gscore < neighbor->gscore_)
                {
                    neighbor->comefrom_ = current;
                    neighbor->gscore_ = tentative_gscore;
                    neighbor->fscore_ = neighbor->gscore_ + GetHeuristic(neighbor);
                    if (!neighbor->traversed_)
                    {
                        neighbor->traversed_ = true;
                        openset_.push(neighbor);
                    }
                }
            }
            
        }
    }

    return false;
}

bool Astar::ReconstructPath(std::vector<Eigen::Vector3d>& path, GridNodePtr current)
{
    
    path.emplace_back(current->coord_);
    GridNodePtr node = current;
    while (node->comefrom_)
    {
        node = node->comefrom_;
        path.emplace_back(node->coord_);
    }
    std::reverse(path.begin(), path.end());
    return true;
}

void Astar::GetNeighborNodes(std::vector<GridNodePtr>& neighbor_nodes, GridNodePtr current)
{
    Eigen::Vector3i idx = current->index_;

    for (auto dir: CONNECTION)
    {
        Eigen::Vector3i n = idx + dir;
        if (!Occupied(n(0), n(1), n(2)))
        {
            neighbor_nodes.emplace_back(grid_map_[n(0)][n(1)][n(2)]);
        }
    }
}

double Astar::GetHeuristic(GridNodePtr current)
{
    current->hscore_ = (current->coord_ - goal_->coord_).norm();
    return current->hscore_;
}

bool Astar::Occupied(int x, int y, int z)
{
    if (x >= 0 && x < GRID_XSIZE && 
        y >= 0 && y < GRID_YSIZE && 
        z >= 0 && z < GRID_ZSIZE)
    {
        return Occupied(grid_map_[x][y][z]);
    }
    else
    {
        return true;
    }
}

inline bool Astar::Occupied(GridNodePtr ptr)
{
    return ptr->occupied_;
}

} /* namespace */