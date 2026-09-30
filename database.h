#pragma once

#include <vector>
#include <string>
#include <map>

struct Element
{
    std::string name;
    double n;
    double A;
};

// ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ (объявление)
extern std::map<std::string, std::vector<double>> E_data;
extern std::map<std::string, std::vector<double>> S_data;
extern std::vector<Element> mixture;

// функции
double interpolate(const std::vector<double>& E, const std::vector<double>& S, double x);
void load_XCOM();
void init_mixture();
double sigma_SM_mix(double E);