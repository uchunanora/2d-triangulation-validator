#ifndef DISPLAY_HPP
#define DISPLAY_HPP

#include <iostream>

template <class T>
void disp(std::string var_name, T &var)
{
    std::cout << var_name << ": " << var << std::endl;
}
template <class T>
void disp_v(vector<T> &vec)
{
    // cout << std::setprecision(16);
    std::cout << "{ ";
    for (const auto &item : vec)
    {
        std::cout << item << " ";
    }
    std::cout << " }" << std::endl;
}
template <class T>
void disp_vv(vector<vector<T>> &vv)
{
    // cout << std::setprecision(16);
    std::cout << "{ \n";
    for (const auto &v : vv)
    {
        for (const auto &item : v)
        {
            std::cout << item << " ";
        }
        std::cout << "\n";
    }
    std::cout << " }" << std::endl;
}

#endif