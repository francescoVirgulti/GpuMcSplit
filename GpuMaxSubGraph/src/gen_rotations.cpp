#include <iostream>
#include <vector>
#include <string>
/*
# generates all possible rotations of a string. Used to find all axes of symmetry in a ring.
# Returns a list with all rotations and respective rotation amounts
*/

std::vector<std::pair<std::string, int> > gen_rotations(const std::string& s);

std::vector<std::pair<std::string, int> > gen_rotations(const std::string& s) {
    int rot_len = s.length();
    std::string tmp = s + s; // Concatenate s with itself
    std::vector<std::pair<std::string, int> > rotations;

    for (int i = 0; i < rot_len; ++i) {
        rotations.push_back(std::make_pair(tmp.substr(i, rot_len), i));
    }
    return rotations;
}