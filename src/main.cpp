// Copyright
#include <cstdarg>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "KDTree.hpp"

#define TEST_BASIC_KD_TREE_ENABLED 1
#define TEST_MODERATE_KD_TREE_ENABLED 1
#define TEST_HARDER_KD_TREE_ENABLED 1
#define TEST_EDGE_CASE_KD_TREE_ENABLED 1
#define TEST_MUTATING_KD_TREE_ENABLED 1
#define TEST_THROWING_KD_TREE_ENABLED 1
#define TEST_CONST_KD_TREE_ENABLED 1


#define TEST_NEAREST_NEIGHBOR_ENABLED 1
#define TEST_MORE_NEAREST_NEIGHBOR_ENABLED 1

#define TEST_BASIC_COPY_ENABLED 1
#define TEST_MODERATE_COPY_ENABLED 1

template<size_t N, typename IteratorType>
Point<N> point_from_range(IteratorType begin, IteratorType end) {
    Point<N> result;
    std::copy(begin, end, result.begin());
    return result;
}

Point<1> make_point(double x) {
    Point<1> result;
    result[0] = x;
    return result;
}

Point<2> make_point(double x, double y) {
    Point<2> result;
    result[0] = x;
    result[1] = y;
    return result;
}

Point<3> make_point(double x, double y, double z) {
    Point<3> result;
    result[0] = x;
    result[1] = y;
    result[2] = z;
    return result;
}

Point<4> make_point(double x, double y, double z, double w) {
    Point<4> result;
    result[0] = x;
    result[1] = y;
    result[2] = z;
    result[3] = w;
    return result;
}

void press_enter_to_continue() {
    std::string line;
    getline(std::cin, line);
}

void do_check_condition(bool expr, const std::string &rationale,
                        const std::string &file, int line) {
    if (expr) {
        std::cout << "PASS: " << rationale << std::endl;
        return;
    }

    std::cout << "FAIL: " << rationale << std::endl;
    std::cout << "  Error at " << file << ", line " << line << std::endl;
    std::cout << "  (ENTER to continue)" << std::endl;

    press_enter_to_continue();
}

void fail_test(const std::exception &e) {
    std::cerr << "TEST FAILED: Unexpected exception: " << e.what() << std::endl;
    press_enter_to_continue();
}

#define CHECK_CONDITION(expr, rationale) \
  do_check_condition(expr, rationale, __FILE__, __LINE__)

void print_banner(const std::string &header) {
    std::cout << "\nBeginning test: " << header << std::endl;
    std::cout << std::setw(40) << std::setfill('-') << "" << std::setfill(' ')
              << std::endl;
}

void test_disabled(const std::string &header) {
    std::cout << "== Test " << header
              << " NOT RUN: press ENTER to continue ==" << std::endl;

    press_enter_to_continue();
}

void end_test() {
    std::cout << "== end of test: press ENTER to continue ==" << std::endl;
    press_enter_to_continue();
}

void test_basic_kd_tree() try {
#if TEST_BASIC_KD_TREE_ENABLED
    print_banner("Basic KDTree Test");

    KDTree<3, size_t> kd;
    CHECK_CONDITION(true, "KDTree construction completed.");

    CHECK_CONDITION(kd.dimension() == 3, "Dimension is three.");
    CHECK_CONDITION(kd.size() == 0, "New KD tree has no elements.");
    CHECK_CONDITION(kd.empty(), "New KD tree is empty.");

    const double data_points[3][3] = {{1, 0, 0},
                                      {0, 1, 0},
                                      {0, 0, 1}};
    for (size_t i = 0; i < 3; ++i) {
        kd.insert(point_from_range<3>(data_points[i], data_points[i] + 3), i);
    }

    CHECK_CONDITION(kd.size() == 3,
                    "After adding three elements, KDTree has size 3.");
    CHECK_CONDITION(!kd.empty(),
                    "After adding three elements, KDTree is not empty.");

    CHECK_CONDITION(
            kd.contains(point_from_range<3>(data_points[0], data_points[0] + 3)),
            "New KD tree has element zero.");
    CHECK_CONDITION(
            kd.contains(point_from_range<3>(data_points[1], data_points[1] + 3)),
            "New KD tree has element one.");
    CHECK_CONDITION(
            kd.contains(point_from_range<3>(data_points[2], data_points[2] + 3)),
            "New KD tree has element two.");


    for (size_t i = 0; i < 3; ++i) {
        CHECK_CONDITION(
                kd.at(point_from_range<3>(data_points[i], data_points[i] + 3)) == i,
                "New KD tree has correct values.");
    }

    end_test();
#else
    test_disabled("test_basic_kd_tree");
#endif
} catch (const std::exception &e) {
    fail_test(e);
}

void test_moderate_kd_tree() try {
#if TEST_MODERATE_KD_TREE_ENABLED
    print_banner("Moderate KDTree Test");

    const double data_points[][4] = {
        {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 0, 1, 1},
        {0, 1, 0, 0}, {0, 1, 0, 1}, {0, 1, 1, 0}, {0, 1, 1, 1},
        {1, 0, 0, 0}, {1, 0, 0, 1}, {1, 0, 1, 0}, {1, 0, 1, 1},
        {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 1, 1, 0}, {1, 1, 1, 1},
    };

    KDTree<4, size_t> kd;
    for (size_t i = 0; i < 16; ++i) {
      kd.insert(point_from_range<4>(data_points[i], data_points[i] + 4), i);
    }

    CHECK_CONDITION(kd.dimension() == 4, "Dimension is four.");
    CHECK_CONDITION(kd.size() == 16,
                    "New KD tree has the right number of elements.");
    CHECK_CONDITION(!kd.empty(), "New KD tree is nonempty.");
    std::cout << "-\n";

    for (size_t i = 0; i < 16; ++i) {
      CHECK_CONDITION(
          kd.at(point_from_range<4>(data_points[i], data_points[i] + 4)) == i,"New KD tree has correct values.");
    }

    std::cout << "-\n";
    CHECK_CONDITION(!kd.contains(make_point(1.0, 1.0, 1.0, 0.5)),
                    "Nonexistent elements aren't in the tree.");
    CHECK_CONDITION(!kd.contains(make_point(0.0, 0.0, 0.0, -0.5)),
                    "Nonexistent elements aren't in the tree.");

    end_test();
#else
    test_disabled("test_moderate_kd_tree");
#endif
} catch (const std::exception &e) {
    fail_test(e);
}

void test_harder_kd_tree() try {
#if TEST_HARDER_KD_TREE_ENABLED
    print_banner("Harder KDTree Test");
    using std::cout, std::endl;

    const double data_points[][4] = {
        {0, 0, 0, 0}, {0, 1, 0, 1}, {0, 0, 0, 0},
        {0, 1, 0, 1}, {0, 1, 1, 0}, {1, 0, 1, 0},
    };
    KDTree<4, size_t> kd;
    cout << "-" << endl;
    for (size_t i = 0; i < 6; ++i) {
      kd.insert(point_from_range<4>(data_points[i], data_points[i] + 4), i);
    }


    CHECK_CONDITION(kd.dimension() == 4, "Dimension is four.");
    CHECK_CONDITION(
        kd.size() == 4,
        "New KD tree has the right number of elements (no duplicates).");
    CHECK_CONDITION(!kd.empty(), "New KD tree is nonempty.");

    for (size_t i = 0; i < 6; ++i) {
      CHECK_CONDITION(
          kd.contains(point_from_range<4>(data_points[i], data_points[i] + 4)),
          "New KD tree has original elems.");
    }

    for (size_t i = 2; i < 6; ++i) {
      CHECK_CONDITION(
          kd.at(point_from_range<4>(data_points[i], data_points[i] + 4)) == i,
          "KD tree has correct labels.");
    }
    for (size_t i = 0; i < 2; ++i) {
      CHECK_CONDITION(
          kd.at(point_from_range<4>(data_points[i], data_points[i] + 4)) == i + 2,
          "Insert overwrites old labels.");
    }

    end_test();
#else
    test_disabled("test_harder_kd_tree");
#endif
} catch (const std::exception &e) {
    fail_test(e);
}

void test_edge_case_kd_tree() try {
#if TEST_EDGE_CASE_KD_TREE_ENABLED
    print_banner("Edge Case KDTree Test");

    const double data_points[8][3] = {
        {0, 0, 0}, {0, 1, 0}, {0, 2, 0}, {0, 3, 0},
        {0, 4, 0}, {0, 5, 0}, {0, 6, 0}, {0, 7, 0},
    };

    KDTree<3, size_t> kd;
    for (size_t i = 0; i < 8; ++i) {
      kd.insert(point_from_range<3>(data_points[i], data_points[i] + 3), i);
    }
    CHECK_CONDITION(kd.dimension() == 3, "Dimension is three.");
    CHECK_CONDITION(kd.size() == 8,
                    "New KD tree has the right number of elements.");
    CHECK_CONDITION(!kd.empty(), "New KD tree is nonempty.");

    for (size_t i = 0; i < kd.size(); ++i) {
      CHECK_CONDITION(
          kd.contains(point_from_range<3>(data_points[i], data_points[i] + 3)),
          "Lookup succeeded.");
    }
    end_test();
#else
    test_disabled("test_edge_case_tree");
#endif
} catch (const std::exception &e) {
    fail_test(e);
}

void test_mutating_kd_tree() try {
#if TEST_MUTATING_KD_TREE_ENABLED
    print_banner("Mutating KDTree Test");

    const double data_points[8][3] = {
        {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
        {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1},
    };

    KDTree<3, size_t> kd;
    for (size_t i = 0; i < 8; ++i) {
      kd[point_from_range<3>(data_points[i], data_points[i] + 3)] = i;
    }
    CHECK_CONDITION(kd.dimension() == 3, "Dimension is three.");
    CHECK_CONDITION(kd.size() == 8,
                    "New KD tree has the right number of elements.");
    CHECK_CONDITION(!kd.empty(), "New KD tree is nonempty.");

    for (size_t i = 0; i < kd.size(); ++i) {
      CHECK_CONDITION(
          kd.contains(point_from_range<3>(data_points[i], data_points[i] + 3)),
          "Lookup succeeded.");
    }
    for (size_t i = 0; i < 8; i += 2) {
      kd[point_from_range<3>(data_points[i], data_points[i] + 3)] = 0;
    }
    for (size_t i = 1; i < 8; i += 2) {
      CHECK_CONDITION(
          kd[point_from_range<3>(data_points[i], data_points[i] + 3)] == i,
          "Keys are correct for odd elements.");
    }
    for (size_t i = 0; i < 8; i += 2) {
      CHECK_CONDITION(
          kd[point_from_range<3>(data_points[i], data_points[i] + 3)] == 0,
          "Keys are correct for even elements.");
    }
    end_test();
#else
    test_disabled("test_mutating_kd_tree");
#endif
} catch (const std::exception &e) {
    fail_test(e);
}

void test_throwing_kd_tree() try {
#if TEST_THROWING_KD_TREE_ENABLED
    print_banner("Throwing KDTree Test");

    {
      KDTree<3, size_t> empty;

      bool did_throw = false;
      try {
        empty.at(make_point(0, 0, 0));
      } catch (const std::out_of_range&) {
        did_throw = true;
      }

      CHECK_CONDITION(did_throw, "Exception generated during non-const lookup.");
    }

    {
      const KDTree<3, size_t> empty;

      bool did_throw = false;
      try {
        empty.at(make_point(0, 0, 0));
      } catch (const std::out_of_range&) {
        did_throw = true;
      }

      CHECK_CONDITION(did_throw, "Exception generated during const lookup.");
    }

    end_test();
#else
    test_disabled("test_throwing_kd_tree");
#endif
} catch (const std::exception &e) {
    fail_test(e);
}

void test_const_kd_tree() try {
#if TEST_CONST_KD_TREE_ENABLED
    print_banner("Const KDTree Test");

    const double data_points[][4] = {
        {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 0, 1, 1}};

    KDTree<4, size_t> kd;
    for (size_t i = 0; i < 4; ++i) {
      kd.insert(point_from_range<4>(data_points[i], data_points[i] + 4), i);
    }

    kd.dimension();
    kd.size();
    kd.empty();
    kd.at(point_from_range<4>(data_points[0], data_points[0] + 4)) = 100;

    const KDTree<4, size_t>& const_kd = kd;

    const_kd.dimension();
    const_kd.size();
    const_kd.empty();
    const_kd.at(point_from_range<4>(data_points[0], data_points[0] + 4));

    CHECK_CONDITION(true, "Const code compiles.");

    CHECK_CONDITION(const_kd.contains(
                        point_from_range<4>(data_points[0], data_points[0] + 4)),
                    "Const KD tree has element zero.");
    CHECK_CONDITION(const_kd.contains(
                        point_from_range<4>(data_points[1], data_points[1] + 4)),
                    "Const KD tree has element one.");
    CHECK_CONDITION(const_kd.contains(
                        point_from_range<4>(data_points[2], data_points[2] + 4)),
                    "Const KD tree has element two.");
    CHECK_CONDITION(const_kd.contains(
                        point_from_range<4>(data_points[3], data_points[3] + 4)),
                    "Const KD tree has element three.");

    CHECK_CONDITION(const_kd.at(point_from_range<4>(data_points[0],
                                                    data_points[0] + 4)) == 100,
                    "Const KD tree has correct values.");
    for (size_t i = 1; i < 4; ++i)
      CHECK_CONDITION(const_kd.at(point_from_range<4>(data_points[i],
                                                      data_points[i] + 4)) == i,
                      "Const KD tree has correct values.");

    end_test();
#else
    test_disabled("test_const_kd_tree");
#endif
} catch (const std::exception &e) {
    std::cout << "Note: KD tree lookup failed, but const code compiles."
              << std::endl;
    fail_test(e);
}

void test_nearest_neighbor() try {
#if TEST_NEAREST_NEIGHBOR_ENABLED
    print_banner("Nearest Neighbor Test");

    const double data_points[][4] = {
        {0, 0, 0, 0},  // 0
        {0, 0, 0, 1},  // 1
        {0, 0, 1, 0},  // 2
        {0, 0, 1, 1},  // 3
        {0, 1, 0, 0},  // 4
        {0, 1, 0, 1},  // 5
        {0, 1, 1, 0},  // 6
        {0, 1, 1, 1},  // 7
        {1, 0, 0, 0},  // 8
        {1, 0, 0, 1},  // 9
        {1, 0, 1, 0},  // 10
        {1, 0, 1, 1},  // 11
        {1, 1, 0, 0},  // 12
        {1, 1, 0, 1},  // 13
        {1, 1, 1, 0},  // 14
        {1, 1, 1, 1},  // 15
    };

    const double testPoints[][4] = {
        {0, 0, 0, 0},           // 0
        {0, 0, 0, 0.7},         // 1
        {0, 0, 0.9, 0},         // 2
        {0, 0, 0.6, 0.6},       // 3
        {0, 0.9, 0, 0},         // 4
        {0, 0.8, 0, 0.7},       // 5
        {0, 0.6, 0.7, -0.1},    // 6
        {-0.4, 0.7, 0.7, 0.7},  // 7
        {1, 0, 0, 0},           // 8
        {1, 0, 0, 1},           // 9
        {1, 0, 1, 0},           // 10
        {1, 0, 1, 1},           // 11
        {1, 1, 0, 0},           // 12
        {1, 1, 0, 1},           // 13
        {1, 1, 1, 0},           // 14
        {1, 1, 1, 1},           // 15
    };

    KDTree<4, size_t> kd;
    for (size_t i = 0; i < 16; ++i)
      kd.insert(point_from_range<4>(data_points[i], data_points[i] + 4), i);

    for (size_t i = 0; i < 16; ++i)
      CHECK_CONDITION(
          kd.knn_value(point_from_range<4>(data_points[i], data_points[i] + 4),
                       1) == i,
          "Nearest neighbor of element is that element.");

    for (size_t i = 0; i < 16; ++i)
      CHECK_CONDITION(
          kd.knn_value(point_from_range<4>(testPoints[i], testPoints[i] + 4),
                       1) == i,
          "Test point yielded correct nearest neighbor.");

    end_test();
#else
    test_disabled("test_nearest_neighbor");
#endif
} catch (const std::exception &e) {
    fail_test(e);
}

void test_more_nearest_neighbor() try {
#if TEST_MORE_NEAREST_NEIGHBOR_ENABLED
    print_banner("More Nearest Neighbor Tests");

    std::vector<std::pair<Point<2>, char> > values;
    values.push_back(std::make_pair(make_point(0.0, 0.0), 'a'));
    values.push_back(std::make_pair(make_point(0.0, 0.5), 'b'));
    values.push_back(std::make_pair(make_point(0.0, 1.0), 'a'));
    values.push_back(std::make_pair(make_point(0.5, 0.0), 'b'));
    values.push_back(std::make_pair(make_point(0.5, 0.5), 'a'));
    values.push_back(std::make_pair(make_point(0.5, 1.0), 'b'));
    values.push_back(std::make_pair(make_point(1.0, 0.0), 'b'));
    values.push_back(std::make_pair(make_point(1.0, 0.5), 'b'));
    values.push_back(std::make_pair(make_point(1.0, 1.0), 'b'));

    KDTree<2, char> kd;
    for (size_t i = 0; i < values.size(); ++i) {
      kd[values[i].first] = values[i].second;
    }

    CHECK_CONDITION(kd.knn_value(make_point(0.0, 0.5), 4) == 'a',
                    "Nearest neighbors are correct.");

    CHECK_CONDITION(kd.knn_value(make_point(0.0, 0.5), 9) == 'b',
                    "Nearest neighbors are correct.");

    CHECK_CONDITION(kd.knn_value(make_point(0.5, 0.5), 5) == 'b',
                    "Nearest neighbors are correct.");

    CHECK_CONDITION(kd.knn_value(make_point(0.75, 0.75), 4) == 'b',
                    "Nearest neighbors are correct.");

    CHECK_CONDITION(kd.knn_value(make_point(+10.0, +10.0), 9) == 'b',
                    "Nearest neighbor at distant point is majority element.");

    CHECK_CONDITION(kd.knn_value(make_point(-10.0, -10.0), 9) == 'b',
                    "No problems with negative values.");

    CHECK_CONDITION(kd.knn_value(make_point(-10.0, -10.0), 25) == 'b',
                    "No problems with looking up more neighbors than elements.");

    end_test();
#else
    test_disabled("test_more_nearest_neighbor");
#endif
} catch (const std::exception &e) {
    fail_test(e);
}

void test_basic_copy() try {
#if TEST_BASIC_COPY_ENABLED
    print_banner("Basic Copy Test");

    KDTree<1, size_t> one;
    for (size_t i = 0; i < 10; ++i) one[make_point(2 * i)] = i;

    {
      KDTree<1, size_t> clone = one;

      CHECK_CONDITION(one.size() == clone.size(),
                      "Copy has the same number of elements as the original.");
      CHECK_CONDITION(one.empty() == clone.empty(),
                      "Copy and original agree on emptiness.");
      CHECK_CONDITION(one.dimension() == clone.dimension(),
                      "Copy and original agree on dimension.");

      for (size_t i = 0; i < 10; ++i)
        CHECK_CONDITION(clone.at(make_point(2 * i)) == i,
                        "Element from original present in copy.");

      for (size_t i = 0; i < 10; ++i)
        CHECK_CONDITION(!clone.contains(make_point(2 * i + 1)),
                        "Other elements not present in copy.");
    }

    {
      KDTree<1, size_t> clone;
      clone = one;

      CHECK_CONDITION(one.size() == clone.size(),
                      "Copy has the same number of elements as the original.");
      CHECK_CONDITION(one.empty() == clone.empty(),
                      "Copy and original agree on emptiness.");
      CHECK_CONDITION(one.dimension() == clone.dimension(),
                      "Copy and original agree on dimension.");

      for (size_t i = 0; i < 10; ++i)
        CHECK_CONDITION(clone.at(make_point(2 * i)) == i,
                        "Element from original present in copy.");

      for (size_t i = 0; i < 10; ++i)
        CHECK_CONDITION(!clone.contains(make_point(2 * i + 1)),
                        "Other elements not present in copy.");
    }

    end_test();
#else
    test_disabled("test_basic_copy");
#endif
} catch (const std::exception &e) {
    fail_test(e);
}

void test_moderate_copy() try {
#if TEST_MODERATE_COPY_ENABLED
    print_banner("Moderate Copy Test");

    KDTree<1, size_t> one;
    for (size_t i = 0; i < 10; ++i) one[make_point(2 * i)] = i;
    {
      KDTree<1, size_t> clone = one;

      for (size_t i = 0; i < 10; ++i) clone[make_point(2 * i + 1)] = i;

      CHECK_CONDITION(one.size() == 10, "Adding to clone change original size.");
      for (size_t i = 0; i < 10; ++i)
        CHECK_CONDITION(!one.contains(make_point(2 * i + 1)),
                        "Modifying copy doesn't modify original.");
    }

    CHECK_CONDITION(one.size() == 10, "After dtor, original is still good.");
    for (size_t i = 0; i < 10; ++i) {
      CHECK_CONDITION(!one.contains(make_point(2 * i + 1)),
                      "After dtor, missing elements still missing.");
      CHECK_CONDITION(one[make_point(2 * i)] == i,
                      "After dtor, original elements are still there.");
    }

    {
      KDTree<1, size_t> clone;
      clone = one;

      clone = clone = (clone = clone);
      (clone = one) = clone;
      clone = clone = clone = clone = clone;
    }

    end_test();
#else
    test_disabled("test_moderate_copy");
#endif
} catch (const std::exception &e) {
    fail_test(e);
}

int main() {

    test_basic_kd_tree();   // PASS
    test_moderate_kd_tree();    // PASS
    test_harder_kd_tree();  // PASS
    test_edge_case_kd_tree();   // PASS
    test_mutating_kd_tree();    // PASS
    test_throwing_kd_tree();    // PASS
    test_const_kd_tree();   // PASS
    test_nearest_neighbor();    // PASS
    test_more_nearest_neighbor();   // PASS
    test_basic_copy();      // PASS
    test_moderate_copy();   // PASS



#if (TEST_BASIC_KD_TREE_ENABLED && TEST_MODERATE_KD_TREE_ENABLED && \
     TEST_HARDER_KD_TREE_ENABLED && TEST_EDGE_CASE_KD_TREE_ENABLED && \
     TEST_MUTATING_KD_TREE_ENABLED && TEST_THROWING_KD_TREE_ENABLED && \
     TEST_CONST_KD_TREE_ENABLED && TEST_NEAREST_NEIGHBOR_ENABLED && \
     TEST_MORE_NEAREST_NEIGHBOR_ENABLED && TEST_BASIC_COPY_ENABLED && \
     TEST_MODERATE_COPY_ENABLED)
    std::cout << "All tests completed!  If they passed, you should be good to go!"
              << std::endl
              << std::endl;
#else
    std::cout
            << "Not all tests were run.  Enable the rest of the tests, then run "
               "again."
            << std::endl
            << std::endl;
#endif

    press_enter_to_continue();

}
