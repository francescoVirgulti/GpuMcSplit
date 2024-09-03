#include <vector>
#include "main.hpp"

// label_classes:    list of all current label classes
// map_size:         size of mapping currently being explored
    // selects label form given label classes, based on the maximum number of nodes in either graph with the
    // specified label.
    // if there is no mapping ignore restriction on adjacency, otherwise return only adjacent lable classes since we
    // are looking for connected sub-graphs. Returns none if there are no label classes adjacent to the current mapping

LabelClass *select_label(std::vector<LabelClass*>& label_classes, int map_size) {
    int min_size = 999;
    LabelClass* label = nullptr;

    for ( LabelClass* c_label : label_classes) {
        if (c_label->adj == 1 || map_size == 0) {
            int c_max_size = std::max(c_label->g.size(), c_label->h.size());
            if (c_max_size < min_size) {
                min_size = c_max_size;
                label = c_label;
            }
        }
    }

    return label;
}