#ifndef LYNDON_NEW_FACTORIZATIONS_H
#define LYNDON_NEW_FACTORIZATIONS_H

#include "finger_graph.h"

Lyndon::factorization cfl(const std::string &s);
Lyndon::factorization icfl(const std::string &s);
Lyndon::factorization cfl_icfl(const std::string &s);
Lyndon::factorization d_cfl(const std::string &s);
Lyndon::factorization d_icfl(const std::string &s);
Lyndon::factorization d_cfl_icfl(const std::string &s);

#endif //LYNDON_NEW_FACTORIZATIONS_H
