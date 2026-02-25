#include "CPMVertexPartition.h"

CPMVertexPartition::CPMVertexPartition(Graph* graph,
      vector<size_t> membership, double resolution_parameter) :
        LinearResolutionParameterVertexPartition(graph,
        membership, resolution_parameter),
        pop_lambda1(0.0), pop_lambda2(0.0), pop_lambda3(0.0),
        pop_lambda4(0.0), pop_lambda5(0.0),
        pop_threshold(0.0), pop_min_threshold(0.0)
{ }

CPMVertexPartition::CPMVertexPartition(Graph* graph,
      vector<size_t> membership) :
        LinearResolutionParameterVertexPartition(graph,
        membership),
        pop_lambda1(0.0), pop_lambda2(0.0), pop_lambda3(0.0),
        pop_lambda4(0.0), pop_lambda5(0.0),
        pop_threshold(0.0), pop_min_threshold(0.0)
{ }

CPMVertexPartition::CPMVertexPartition(Graph* graph,
      double resolution_parameter) :
        LinearResolutionParameterVertexPartition(graph, resolution_parameter),
        pop_lambda1(0.0), pop_lambda2(0.0), pop_lambda3(0.0),
        pop_lambda4(0.0), pop_lambda5(0.0),
        pop_threshold(0.0), pop_min_threshold(0.0)
{ }

CPMVertexPartition::CPMVertexPartition(Graph* graph) :
        LinearResolutionParameterVertexPartition(graph),
        pop_lambda1(0.0), pop_lambda2(0.0), pop_lambda3(0.0),
        pop_lambda4(0.0), pop_lambda5(0.0),
        pop_threshold(0.0), pop_min_threshold(0.0)
{ }

CPMVertexPartition::CPMVertexPartition(Graph* graph,
      vector<size_t> membership, double resolution_parameter,
      double pop_lambda1, double pop_lambda2, double pop_lambda3,
      double pop_lambda4, double pop_lambda5,
      double pop_threshold, double pop_min_threshold) :
        LinearResolutionParameterVertexPartition(graph,
        membership, resolution_parameter),
        pop_lambda1(pop_lambda1), pop_lambda2(pop_lambda2), pop_lambda3(pop_lambda3),
        pop_lambda4(pop_lambda4), pop_lambda5(pop_lambda5),
        pop_threshold(pop_threshold), pop_min_threshold(pop_min_threshold)
{ }

CPMVertexPartition::CPMVertexPartition(Graph* graph,
      double resolution_parameter,
      double pop_lambda1, double pop_lambda2, double pop_lambda3,
      double pop_lambda4, double pop_lambda5,
      double pop_threshold, double pop_min_threshold) :
        LinearResolutionParameterVertexPartition(graph, resolution_parameter),
        pop_lambda1(pop_lambda1), pop_lambda2(pop_lambda2), pop_lambda3(pop_lambda3),
        pop_lambda4(pop_lambda4), pop_lambda5(pop_lambda5),
        pop_threshold(pop_threshold), pop_min_threshold(pop_min_threshold)
{ }

CPMVertexPartition::~CPMVertexPartition()
{ }

CPMVertexPartition* CPMVertexPartition::create(Graph* graph)
{
  return new CPMVertexPartition(graph, this->resolution_parameter,
                                this->pop_lambda1, this->pop_lambda2, this->pop_lambda3,
                                this->pop_lambda4, this->pop_lambda5,
                                this->pop_threshold, this->pop_min_threshold);
}

CPMVertexPartition* CPMVertexPartition::create(Graph* graph, vector<size_t> const& membership)
{
  return new CPMVertexPartition(graph, membership, this->resolution_parameter,
                                this->pop_lambda1, this->pop_lambda2, this->pop_lambda3,
                                this->pop_lambda4, this->pop_lambda5,
                                this->pop_threshold, this->pop_min_threshold);
}

/********************************************************************************
  RBER implementation of a vertex partition
  (which includes a resolution parameter).
 ********************************************************************************/
double CPMVertexPartition::diff_move(size_t v, size_t new_comm)
{
  #ifdef DEBUG
    cerr << "double CPMVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
    cerr << "Using resolution parameter: " << this->resolution_parameter << "." << endl;
  #endif
  size_t old_comm = this->membership(v);
  double diff = 0.0;
  if (new_comm != old_comm)
  {
    double w_to_old = this->weight_to_comm(v, old_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_to_old: " << w_to_old << endl;
    #endif
    double w_to_new = this->weight_to_comm(v, new_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_to_new: " << w_to_new << endl;
    #endif
    double w_from_old = this->weight_from_comm(v, old_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_from_old: " << w_from_old << endl;
    #endif
    double w_from_new = this->weight_from_comm(v, new_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_from_new: " << w_from_new << endl;
    #endif
    double nsize = this->graph->node_size(v);
    #ifdef DEBUG
      cerr << "\t" << "nsize: " << nsize << endl;
    #endif
    double csize_old = this->csize(old_comm);
    #ifdef DEBUG
      cerr << "\t" << "csize_old: " << csize_old << endl;
    #endif
    double csize_new = this->csize(new_comm);
    #ifdef DEBUG
      cerr << "\t" << "csize_new: " << csize_new << endl;
    #endif
    double self_weight = this->graph->node_self_weight(v);
    #ifdef DEBUG
      cerr << "\t" << "self_weight: " << self_weight << endl;
      cerr << "\t" << "density: " << this->graph->density() << endl;
    #endif
    double possible_edge_difference_old = 0.0;
    if (this->graph->correct_self_loops())
      possible_edge_difference_old = nsize*(2.0*csize_old - nsize);
    else
      possible_edge_difference_old = nsize*(2.0*csize_old - nsize - 1.0);
    #ifdef DEBUG
      cerr << "\t" << "possible_edge_difference_old: " << possible_edge_difference_old << endl;
    #endif
    double diff_old = w_to_old + w_from_old -
        self_weight - this->resolution_parameter*possible_edge_difference_old;
    #ifdef DEBUG
      cerr << "\t" << "diff_old: " << diff_old << endl;
    #endif
    double possible_edge_difference_new = 0.0;
    if (this->graph->correct_self_loops())
      possible_edge_difference_new = nsize*(2.0*csize_new + nsize);
    else
      possible_edge_difference_new = nsize*(2.0*csize_new + nsize - 1.0);
    #ifdef DEBUG
      cerr << "\t" << "possible_edge_difference_new: " << possible_edge_difference_new << endl;
    #endif
    double diff_new = w_to_new + w_from_new + self_weight -
        this->resolution_parameter*possible_edge_difference_new;
    #ifdef DEBUG
      cerr << "\t" << "diff_new: " << diff_new << endl;
    #endif
    diff = diff_new - diff_old;
    #ifdef DEBUG
      cerr << "\t" << "diff: " << diff << endl;;
    #endif

    // Population penalty
    if (this->pop_lambda1 > 0.0 || this->pop_lambda2 > 0.0 || this->pop_lambda3 > 0.0 ||
        this->pop_lambda4 > 0.0 || this->pop_lambda5 > 0.0)
    {
      double node_pop       = this->graph->node_pop(v);
      double pop_old_before = this->cpop(old_comm);
      double pop_old_after  = pop_old_before - node_pop;
      double pop_new_before = this->cpop(new_comm);
      double pop_new_after  = pop_new_before + node_pop;

      double penalty_before = 0.0;
      double penalty_after  = 0.0;

      // Upper bound: polynomial penalty for communities that exceed pop_threshold
      if (this->pop_threshold > 0.0)
      {
        // Calculate penalty_before for old community
        if (pop_old_before > this->pop_threshold)
        {
          double excess = pop_old_before - this->pop_threshold;
          double term = excess;
          penalty_before += this->pop_lambda1 * term;
          term *= excess;
          penalty_before += this->pop_lambda2 * term;
          term *= excess;
          penalty_before += this->pop_lambda3 * term;
          term *= excess;
          penalty_before += this->pop_lambda4 * term;
          term *= excess;
          penalty_before += this->pop_lambda5 * term;
        }

        // Calculate penalty_before for new community
        if (pop_new_before > this->pop_threshold)
        {
          double excess = pop_new_before - this->pop_threshold;
          double term = excess;
          penalty_before += this->pop_lambda1 * term;
          term *= excess;
          penalty_before += this->pop_lambda2 * term;
          term *= excess;
          penalty_before += this->pop_lambda3 * term;
          term *= excess;
          penalty_before += this->pop_lambda4 * term;
          term *= excess;
          penalty_before += this->pop_lambda5 * term;
        }

        // Calculate penalty_after for old community
        if (pop_old_after > this->pop_threshold)
        {
          double excess = pop_old_after - this->pop_threshold;
          double term = excess;
          penalty_after += this->pop_lambda1 * term;
          term *= excess;
          penalty_after += this->pop_lambda2 * term;
          term *= excess;
          penalty_after += this->pop_lambda3 * term;
          term *= excess;
          penalty_after += this->pop_lambda4 * term;
          term *= excess;
          penalty_after += this->pop_lambda5 * term;
        }

        // Calculate penalty_after for new community
        if (pop_new_after > this->pop_threshold)
        {
          double excess = pop_new_after - this->pop_threshold;
          double term = excess;
          penalty_after += this->pop_lambda1 * term;
          term *= excess;
          penalty_after += this->pop_lambda2 * term;
          term *= excess;
          penalty_after += this->pop_lambda3 * term;
          term *= excess;
          penalty_after += this->pop_lambda4 * term;
          term *= excess;
          penalty_after += this->pop_lambda5 * term;
        }
      }

      // Lower bound: keep existing simple linear penalty for lower threshold
      if (this->pop_min_threshold > 0.0 && this->pop_lambda1 > 0.0)
      {
        // old community before move (always non-empty since v is in it)
        if (pop_old_before < this->pop_min_threshold)
          penalty_before += this->pop_lambda1 * (this->pop_min_threshold - pop_old_before);
        // old community after move (only penalise if it still has members)
        if (pop_old_after > 0.0 && pop_old_after < this->pop_min_threshold)
          penalty_after += this->pop_lambda1 * (this->pop_min_threshold - pop_old_after);

        // new community before move (only penalise if it has members)
        if (pop_new_before > 0.0 && pop_new_before < this->pop_min_threshold)
          penalty_before += this->pop_lambda1 * (this->pop_min_threshold - pop_new_before);
        // new community after move (always has members since v just joined)
        if (pop_new_after < this->pop_min_threshold)
          penalty_after += this->pop_lambda1 * (this->pop_min_threshold - pop_new_after);
      }

      diff -= (penalty_after - penalty_before);
    }
  }
  #ifdef DEBUG
    cerr << "exit CPMVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
    cerr << "return " << diff << endl << endl;
  #endif
  return diff;
}

double CPMVertexPartition::quality(double resolution_parameter)
{
  #ifdef DEBUG
    cerr << "double CPMVertexPartition::quality()" << endl;
  #endif
  double mod = 0.0;
  for (size_t c = 0; c < this->n_communities(); c++)
  {
    double csize = this->csize(c);
    double w = this->total_weight_in_comm(c);
    double comm_possible_edges = this->graph->possible_edges(csize);

    #ifdef DEBUG
      cerr << "\t" << "Comm: " << c << ", w_c=" << w << ", n_c=" << csize << ", comm_possible_edges=" << comm_possible_edges << ", p=" << this->graph->density() << "." << endl;
    #endif
    mod += w - resolution_parameter*comm_possible_edges;

    // Population penalty
    if (this->pop_lambda1 > 0.0 || this->pop_lambda2 > 0.0 || this->pop_lambda3 > 0.0 ||
        this->pop_lambda4 > 0.0 || this->pop_lambda5 > 0.0)
    {
      double cpop = this->cpop(c);
      // Upper bound: polynomial penalty
      if (this->pop_threshold > 0.0 && cpop > this->pop_threshold)
      {
        double excess = cpop - this->pop_threshold;
        double term = excess;
        mod -= this->pop_lambda1 * term;
        term *= excess;
        mod -= this->pop_lambda2 * term;
        term *= excess;
        mod -= this->pop_lambda3 * term;
        term *= excess;
        mod -= this->pop_lambda4 * term;
        term *= excess;
        mod -= this->pop_lambda5 * term;
      }
      // Lower bound (only for non-empty communities): simple linear penalty
      if (this->pop_min_threshold > 0.0 && cpop > 0.0 && cpop < this->pop_min_threshold && this->pop_lambda1 > 0.0)
        mod -= this->pop_lambda1 * (this->pop_min_threshold - cpop);
    }
  }
  #ifdef DEBUG
    cerr << "exit double CPMVertexPartition::quality()" << endl;
    cerr << "return " << mod << endl << endl;
  #endif
  return (2.0 - this->graph->is_directed())*mod;
}

