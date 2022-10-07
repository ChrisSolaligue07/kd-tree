// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"
#include "BoundedPQueue.h"

template<size_t N, typename ElemType>
class KDTree {
public:
    struct KdNode {
        // ElemType -> Tipo de valor del dato del nodo.
        // N -> Cantidad de dimensiones por punto.
        Point<N> punto_s;
        KdNode *moves[2];
        ElemType value;
        size_t nivel = 0;

        KdNode(Point<N> punto_s, ElemType value) {
            moves[0] = moves[1] = nullptr;
            this->punto_s = punto_s;
            this->value = value;
        }
    };

    typedef std::pair<Point<N>, ElemType> value_type;

    KDTree();

    ~KDTree();

    KDTree(const KDTree &rhs);

    KDTree &operator=(const KDTree &rhs);

    [[nodiscard]] size_t dimension() const;

    [[nodiscard]] size_t size() const;

    [[nodiscard]] bool empty() const;

    bool contains(const Point<N> &pt) const;

    bool find(Point<N> pt, KdNode *&node_ant, KdNode *&node_sig, int &axis);

    bool find(KdNode *&i, Point<N> pt) const;

    bool insert(const Point<N> &pt, const ElemType &value);

    ElemType &operator[](const Point<N> &pt);

    ElemType &at(const Point<N> &pt);

    const ElemType &at(const Point<N> &pt) const;

    ElemType knn_value(const Point<N> &key, size_t k) const;

    std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;

    void nearest_neighbor(KdNode *current_node, const Point<N> &key, BoundedPQueue<ElemType> &nearest_neighbors_candidates, int depth) const;

    void _nearest_neighbor_(KdNode *&current_node, KdNode *&nearest_neighbors_candidate, const Point<N> &key, double &best_distance, int depth) const;

    KdNode *cpy_tree(KdNode *_root_, Point<N> &pt, const ElemType &value);

    void d_resource(KdNode *current_node);

    KdNode *root;
private:

    size_t dimension_;
    size_t size_;

};

template<size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    root = nullptr;
    dimension_ = N;
    size_ = 0;
}

template<size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    KdNode *tmp;
    for(KdNode *i = root; i; i = tmp){
        if(!i->moves[0]){
            tmp = i->moves[1];
            delete i;
        } else {
            tmp = i->moves[0];
            i->moves[0] = tmp->moves[1];
            tmp->moves[1] = i;
        }
    }
}

template<size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
    // TODO(me): Fill this in.
    root = cpy_tree(rhs.root, rhs.root->punto_s, rhs.root->value);
    dimension_ = N;
    size_ = rhs.size_;
}

template<size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
    // TODO(me): Fill this in.
    if (this != &rhs){
        d_resource(root);
        root = cpy_tree(rhs.root, rhs.root->punto_s, rhs.root->value);
        size_ = rhs.size_;
    }
    return *this;
}

template<size_t N, typename ElemType>
typename KDTree<N, ElemType>::KdNode *KDTree<N, ElemType>::cpy_tree(KdNode *_root_, Point<N> &pt, const ElemType &value) {
    if (_root_ == NULL){
        return NULL;
    }
    KdNode *n_node = new KdNode(_root_->punto_s, _root_->value);
    n_node->moves[0] = cpy_tree(_root_->moves[0], _root_->punto_s, _root_->value);
    n_node->moves[1] = cpy_tree(_root_->moves[1], _root_->punto_s, _root_->value);
    return n_node;
}

template<size_t N, typename ElemType>
void KDTree<N, ElemType>::d_resource(KdNode *current_node) {
    if (!current_node) {
        return;
    }
    d_resource(current_node->moves[0]);
    d_resource(current_node->moves[1]);
    delete current_node;
}
// ------------------------------------------------------------------------------

template<size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
    return N;
}

template<size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    return size_;
}

template<size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
    if (size_ == 0){
        return true;
    } else {
        return false;
    }
}

template<size_t N, typename ElemType>
bool KDTree<N, ElemType>::find(KdNode *&i, Point<N> pt) const {
    int axis = 0;
    for(i = root ; i && i->punto_s != pt ;i = i->moves[pt[axis] >= i->punto_s[axis]], axis++){
        if (axis == pt.size()) { axis = 0; }
    }
    if (i != nullptr && i->punto_s == pt) {
        return true;
    } else {
        return false;
    }
}

template<size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
    /*
    int axis = 0;
    for(KdNode *i = root; i; i = i->moves[ pt[axis] >= i->punto_s[axis] ], axis++){
        if (pt == i->punto_s){
            return true;
        }
        if (axis == pt.size()){ axis = 0; }
    }
    return false;
    */

    /*if (find(pt)){
        return true;
    } else {
        return false;
    }*/
    KdNode *i;
    return find(i, pt);
}

template<size_t N, typename ElemType>
bool KDTree<N, ElemType>::find(const Point<N> pt, KDTree::KdNode *&node_ant, KDTree::KdNode *&node_sig, int &axis){
    // Una dimensión -> Inserción Lineal
    // Más de dos dimensiones -> Inserción por axis
    // pt : dato de que se busca
    // root->punto_s[] : dato dentro del árbol

    for(node_sig = root; node_sig && node_sig->punto_s != pt
            ; node_sig = node_sig->moves[pt[axis] >= node_sig->punto_s[axis]], axis++){
        if(axis == pt.size()){ axis = 0; }
        node_ant = node_sig;
    }
    if(node_sig != nullptr && node_sig->punto_s == pt){
        // std::cout << "node_sig: " << node_sig->value << std::endl;
        return false;
    }
    if (node_ant != nullptr && node_ant->punto_s != pt) {
        // std::cout << "node_ant: " << node_ant->value << std::endl;
        return true;
    } else {
        return false;
    }
}

template<size_t N, typename ElemType>
bool KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
    if (root == nullptr) {
        root = new KdNode(pt, value);
        size_ = 1;
        return true;
    } else {
        KdNode *node_ant, *node_sig;
        int axis = 0;
        if (find(pt, node_ant, node_sig, axis)) {
            axis--;
            // std::cout << "Axis : " << axis << std::endl;
            node_sig = new KdNode(pt, value);
            node_sig->nivel ++;
            if (pt[axis] >= node_ant->punto_s[axis]) {
                node_ant->moves[1] = node_sig;
                size_++;
            } else {
                node_ant->moves[0] = node_sig;
                size_++;
            }
            return true;
        } else {
            if(node_sig != nullptr){
                // std::cout << "Replace data..." << std::endl;
                node_sig->value = value;    // Valor ya existe, solo se actualiza su valor
            }
            return false;
        }
    }
}

template<size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
    /*
     *  Point<N>& pt);	Retorna una referencia al valor asociado con el punto pt.
     *  Si el punto no existe en el KD-Tree, su valor es añadido con el valor por defecto de ElemenType
     *  , y una referencia a su nuevo valor es retornado.
     *  Este es el mismo comportamiento del map operator[] STL.
     *  Tengan en cuenta que esta función no sobrecarga const por que el árbol puede mutar.
     */
    KdNode *i = root;
    /*int axis = 0;
    for(; i; i = i->moves[pt[axis] >= i->punto_s[axis]], axis++){
        if (i->punto_s == pt){
            return (i->value);  // Si el punto existe retorna el valor
        }
        if (axis == pt.size()){ axis = 0; }
    }
    // std::cout << "Punto no existe - Creando uno nuevo." << std::endl;
    if (i == nullptr){
        insert(pt, {});
        i = root;
        axis = 0;
        for(; i; i = i->moves[pt[axis] >= i->punto_s[axis]], axis++){
            if (i->punto_s == pt){
                return (i->value);
            }
            if (axis == pt.size()){ axis = 0; }
        }
    }
*/
    if (find(i, pt)){
        return i->value;
    } else {
        // std::cout << "Punto no existe - Creando uno nuevo." << std::endl;
        insert(pt, {});
        find(i, pt);
        return i->value;
    }
}

template<size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
    const auto &th = *this;
    return const_cast<ElemType&>(th.at(pt));
}

template<size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
    int axis = 0;
    for(KdNode *i = root; i; i = i->moves[pt[axis] >= i->punto_s[axis]], axis++){
        if (pt == i->punto_s){
            return i->value;
        }
        if (axis == pt.size()){ axis = 0; }
    }
    throw std::out_of_range("No se encontro el punto.");
}

template<size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
    // TODO(me): Fill this in.
    /*
     * Dado un punto v y un entero key, encuentra los k puntos en el KD-Tree más cercanos a key y
     * devuelve el valor más común asociado con esos puntos. En caso de empate, se elegirá uno de
     * los valores más frecuentes.
     */
    BoundedPQueue<ElemType> pQueue(k);
    if (empty())
        return ElemType();

    nearest_neighbor(root, key, pQueue, {});

/*
    KdNode *n = nullptr;
    KdNode *curr = root;
    double best = std::numeric_limits<double>::infinity();
    _nearest_neighbor_(curr, n, key, best, {});
*/

    std::unordered_map<ElemType, int> counter;
    while (!pQueue.empty()) {
        ++counter[pQueue.dequeueMin()];
    }

    ElemType result;
    int cnt = -1;
    for (const auto &p : counter) {
        if (p.second > cnt) {
            result = p.first;
            cnt = p.second;
        }
    }
/*
    std::cout << "A : " << result << "\n";
    std::cout << "B : " << n->value << "\n";
*/

    return result;
}

template<size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key, size_t k) const {
    // TODO(me): Fill this in.
    // Retorna los k valores más cercanos al punto key contenidos en un vector.
    std::vector<ElemType> values;
    return values;
}

template<size_t N, typename ElemType>
void KDTree<N, ElemType>::nearest_neighbor(KdNode *current_node, const Point<N> &key, BoundedPQueue<ElemType> &nearest_neighbors_candidates, int depth) const {
    if (!current_node){
        return;
    }
    nearest_neighbors_candidates.enqueue(current_node->value, (distance(current_node->punto_s, key)));
    int axis = depth % dimension_;
    bool right = false;
    if (key[axis] < current_node->punto_s[axis]){
        right = true;
        nearest_neighbor(current_node->moves[0], key, nearest_neighbors_candidates, ++depth);
    } else {
        right = false;
        nearest_neighbor(current_node->moves[1], key, nearest_neighbors_candidates, ++depth);
    }
    if (nearest_neighbors_candidates.size() < nearest_neighbors_candidates.maxSize() || fabs(current_node->punto_s[axis] - key[axis]) < nearest_neighbors_candidates.worst()){
        if (right) {
            nearest_neighbor(current_node->moves[1], key, nearest_neighbors_candidates, ++depth);
        } else {
            nearest_neighbor(current_node->moves[0], key, nearest_neighbors_candidates, ++depth);
        }
    }
}

template<size_t N, typename ElemType>
void KDTree<N, ElemType>::_nearest_neighbor_(KdNode *&current_node, KdNode *&nearest_neighbors_candidate,
                                             const Point<N> &key, double &best_distance, int depth) const {
    if (!current_node){
        return;
    }
    if (distance(current_node->punto_s, key) < best_distance){
        best_distance = distance(current_node->punto_s, key);
        nearest_neighbors_candidate = current_node;
    }
    int axis = depth % dimension_;
    bool right = false;
    if (key[axis] < current_node->punto_s[axis]){
        right = true;
        _nearest_neighbor_(current_node->moves[0], nearest_neighbors_candidate, key, best_distance, ++depth);
    } else {
        right = false;
        _nearest_neighbor_(current_node->moves[1], nearest_neighbors_candidate, key, best_distance, ++depth);
    }

    if (fabs(current_node->punto_s[axis] - key[axis]) < distance(nearest_neighbors_candidate->punto_s, key)){
        // if (fabs(current_node->punto_s[axis] - key[axis]) < std::numeric_limits<double>::infinity()){
        if (right){
            _nearest_neighbor_(current_node->moves[1], nearest_neighbors_candidate, key, best_distance, ++depth);
        } else {
            _nearest_neighbor_(current_node->moves[0], nearest_neighbors_candidate, key, best_distance, ++depth);
        }
    }
}

#endif  // SRC_KDTREE_HPP_
