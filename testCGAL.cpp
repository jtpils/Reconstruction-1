// Constructing an arrangement using the simple edge-insertion functions.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include "arr_print.h"
typedef int                                           Number_type;
typedef CGAL::Simple_cartesian<Number_type>           Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;
int main()
{
    Arrangement_2   arr;
    Segment_2 s1(Point_2(0, 3), Point_2(10, 3));
    insert(arr, s1);
    Segment_2 s2(Point_2(10, 7), Point_2(0, 7));
    insert(arr, s2);
    Segment_2 s22(Point_2(0, 8), Point_2(10, 8));
    insert(arr, s22);
    Segment_2 s3(Point_2(3, 0), Point_2(3, 10));
    insert(arr, s3);
    Segment_2 s4(Point_2(7, 0), Point_2(7, 10));
    insert(arr, s4);

    // Print the arrangement faces.
    Arrangement_2::Face_const_iterator    fit;
    std::cout << arr.number_of_faces() << " faces:" << std::endl;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
        // Print the outer boundary.
        if (fit->is_unbounded())
            continue;
        else
        {
            std::cout << "Outer boundary: ";
            print_ccb<Arrangement_2> (fit->outer_ccb());
        }
        // Print the boundary of each of the holes.
        Arrangement_2::Hole_const_iterator  hole;
        int                                         index = 1;
        for (hole = fit->holes_begin(); hole != fit->holes_end(); ++hole, ++index)
        {
            std::cout << "    Hole #" << index << ": ";
            print_ccb<Arrangement_2> (*hole);
        }
        // Print the isolated vertices.
        Arrangement_2::Isolated_vertex_const_iterator  iv;
        for (iv = fit->isolated_vertices_begin(), index = 1;
             iv != fit->isolated_vertices_end(); ++iv, ++index)
        {
            std::cout << "    Isolated vertex #" << index << ": "
                      << "(" << iv->point() << ")" << std::endl;
        }
    return 0;
}