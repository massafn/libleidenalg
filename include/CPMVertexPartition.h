#ifndef CPMVERTEXPARTITION_H
#define CPMVERTEXPARTITION_H

#include "LinearResolutionParameterVertexPartition.h"

class LIBLEIDENALG_EXPORT CPMVertexPartition : public LinearResolutionParameterVertexPartition
{
  public:
    CPMVertexPartition(Graph* graph,
          vector<size_t> membership, double resolution_parameter);
    CPMVertexPartition(Graph* graph,
          vector<size_t> membership);
    CPMVertexPartition(Graph* graph,
      double resolution_parameter);
    CPMVertexPartition(Graph* graph);

    // Constructors with population penalty parameters
    CPMVertexPartition(Graph* graph,
          vector<size_t> membership, double resolution_parameter,
          double pop_lambda, double pop_threshold,
          double pop_min_threshold = 0.0);
    CPMVertexPartition(Graph* graph,
          double resolution_parameter,
          double pop_lambda, double pop_threshold,
          double pop_min_threshold = 0.0);

    virtual ~CPMVertexPartition();
    virtual CPMVertexPartition* create(Graph* graph);
    virtual CPMVertexPartition* create(Graph* graph, vector<size_t> const& membership);

    virtual double diff_move(size_t v, size_t new_comm);
    virtual double quality(double resolution_parameter);

    double pop_lambda;         // penalty weight
    double pop_threshold;      // population upper threshold (0 = disabled)
    double pop_min_threshold;  // population lower threshold (0 = disabled)

  protected:
  private:
};

#endif // CPMVERTEXPARTITION_H
