//
// Created by Henry Hollis on 1/27/24.
//

#ifndef LEIDEN_CCD_CCDMODULARITYVERTEXPARTITION_H
#define LEIDEN_CCD_CCDMODULARITYVERTEXPARTITION_H
#define CCD_COMM_SIZE 2
#include "MutableVertexPartition.h"
typedef vector<double> Vector;

class LIBLEIDENALG_EXPORT ccdModularityVertexPartition : public MutableVertexPartition
{
public:
    //constructors
    // Constructor with matrix initialization
    ccdModularityVertexPartition(Graph* graph,
                                 vector<size_t> const& membership,
                                 const std::vector<Vector>& geneSampleMatrix);
    ccdModularityVertexPartition(Graph* graph,
                              vector<size_t> const& membership);
    ccdModularityVertexPartition(Graph* graph);

    virtual ~ccdModularityVertexPartition();
    virtual ccdModularityVertexPartition* create(Graph* graph);
    virtual ccdModularityVertexPartition* create(Graph* graph, vector<size_t> const& membership);
    virtual ccdModularityVertexPartition* create(Graph* graph, vector<size_t> const& membership, const std::vector<Vector>& geneSampleMatrix);

    virtual double diff_move(size_t v, size_t new_comm);
    virtual double quality();
    //  Setter method for the matrix
    void setGeneSampleMatrix(const std::vector<Vector>& geneSampleMatrix);

    // Getter for geneSampleMatrix
    [[nodiscard]] const std::vector<Vector>& getMatrix();

protected:
private:
    // Matrix representing genes and samples
    std::vector<Vector> geneSampleMatrix;

};

#endif //LEIDEN_CCD_CCDMODULARITYVERTEXPARTITION_H
