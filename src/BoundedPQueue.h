//
// Created by chris on 6/10/2022.
//

#ifndef KDTREE_BOUNDEDPQUEUE_H
#define KDTREE_BOUNDEDPQUEUE_H

#include <map>
#include <algorithm>
#include <limits>

template <typename T>
class BoundedPQueue {
public:
    // Constructor: BoundedPQueue(size_t maxSize);
    // Usage: BoundedPQueue<int> bpq(15);
    // --------------------------------------------------
    // Constructs a new, empty BoundedPQueue with
    // maximum size equal to the constructor argument.
    ///
    explicit BoundedPQueue(std::size_t maxSize);

    // void enqueue(const T& value, double priority);
    // Usage: bpq.enqueue("Hi!", 2.71828);
    // --------------------------------------------------
    // Enqueues a new element into the BoundedPQueue with
    // the specified priority. If this overflows the maximum
    // size of the queue, the element with the highest
    // priority will be deleted from the queue. Note that
    // this might be the element that was just added.
    void enqueue(const T& value, double priority);

    // T dequeueMin();
    // Usage: int val = bpq.dequeueMin();
    // --------------------------------------------------
    // Returns the element from the BoundedPQueue with the
    // smallest priority value, then removes that element
    // from the queue.
    T dequeueMin();

    // size_t size() const;
    // bool empty() const;
    // Usage: while (!bpq.empty()) { ... }
    // --------------------------------------------------
    // Returns the number of elements in the queue and whether
    // the queue is empty, respectively.
    std::size_t size() const;
    bool empty() const;

    // size_t maxSize() const;
    // Usage: size_t queueSize = bpq.maxSize();
    // --------------------------------------------------
    // Returns the maximum number of elements that can be
    // stored in the queue.
    std::size_t maxSize() const;

    // double best() const;
    // double worst() const;
    // Usage: double highestPriority = bpq.worst();
    // --------------------------------------------------
    // best() returns the smallest priority of an element
    // stored in the container (i.e. the priority of the
    // element that will be dequeued first using dequeueMin).
    // worst() returns the largest priority of an element
    // stored in the container.  If an element is enqueued
    // with a priority above this value, it will automatically
    // be deleted from the queue.  Both functions return
    // numeric_limits<double>::infinity() if the queue is
    // empty.
    double best()  const;
    double worst() const;

private:
    // This class is layered on top of a multimap mapping from priorities
    // to elements with those priorities.
    std::multimap<double, T> elems;
    std::size_t maximumSize;
};

/** BoundedPQueue class implementation details */

template <typename T>
BoundedPQueue<T>::BoundedPQueue(std::size_t maxSize) {
    maximumSize = maxSize;
}

// enqueue adds the element to the map, then deletes the last element of the
// map if there size exceeds the maximum size.
template <typename T>
void BoundedPQueue<T>::enqueue(const T& value, double priority) {
    // Add the element to the collection.
    elems.insert(std::make_pair(priority, value));

    // If there are too many elements in the queue, drop off the last one.
    if (size() > maxSize()) {
        typename std::multimap<double, T>::iterator last = elems.end();
        --last; // Now points to highest-priority element
        elems.erase(last);
    }
}

// dequeueMin copies the lowest element of the map (the one pointed at by
// begin()) and then removes it.
template <typename T>
T BoundedPQueue<T>::dequeueMin() {
    // Copy the best value.
    T result = elems.begin()->second;

    // Remove it from the map.
    elems.erase(elems.begin());

    return result;
}

// size() and empty() call directly down to the underlying map.
template <typename T>
std::size_t BoundedPQueue<T>::size() const {
    return elems.size();
}

template <typename T>
bool BoundedPQueue<T>::empty() const {
    return elems.empty();
}

// maxSize just returns the appropriate data member.
template <typename T>
std::size_t BoundedPQueue<T>::maxSize() const {
    return maximumSize;
}

// The best() and worst() functions check if the queue is empty,
// and if so return infinity.
template <typename T>
double BoundedPQueue<T>::best() const {
    return empty()? std::numeric_limits<double>::infinity() : elems.begin()->first;
}

template <typename T>
double BoundedPQueue<T>::worst() const {
    return empty()? std::numeric_limits<double>::infinity() : elems.rbegin()->first;
}


#endif //KDTREE_BOUNDEDPQUEUE_H
