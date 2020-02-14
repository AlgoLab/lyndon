#include <fstream>
#include <iostream>
#include "finger_graph.h"
#include "utils.h"
#include "argagg.h"
#include <zlib.h>
#include <ctime>

using namespace std;
using namespace Lyndon;

void print_time() {
    time_t t = time(0);
    tm *ltm = localtime(&t);
    fprintf(stderr, "[%d:%d] - ", ltm->tm_hour, ltm->tm_min);
}

//  N.B. A differenza dell'implementazione in Python, il file in input e' nel formato
//      `read_id` `offset`|`f1` `f2` `f3`...
//  Ai fattori e' stato gia' rimosso il bordo (oltre che applicato l'algoritmo) e `offset` e' il numero di basi rimosso.
//  Stessa cosa vale per il file delle fattorizzazioni.
int main(int argc, char *argv[]) {
    argagg::parser argparser {{
        { "help", {"-h", "--help"},
        "help", 0},
        { "k", {"-k"},
        "k-finger dimension", 1},
        { "limit", {"-l", "--limit"},
        "minimum k-finger length", 1},
        { "no_norm", {"--no-norm"},
        "do not normalize k-fingers", 0},
        { "no_enriched", {"--no-enriched"},
          "do not enrich k-fingers", 0},
    }};

    argagg::parser_results args;
    try {
        args = argparser.parse(argc, argv);
    } catch (const std::exception& e) {
        cerr << e.what() << endl;
        return 1;
    }

    ostringstream usage;
    usage << "Usage: " << argv[0] << " [-k k] [-l limit] [--no-norm] [--no-enriched] FACTORS_PATH" << endl << endl;
    if (args["help"]) {
        cerr << usage.str();
        return 0;
    }

    if (args.pos.size() == 0) {
        cerr << usage.str();
        return 1;
    }

    auto k = args["k"].as<int>(5);
    auto limit = args["limit"].as<int>(30);
    bool no_norm = args["no_norm"];
    bool no_enriched = args["no_enriched"];
    auto factors_path = args.pos[0];

    ifstream in(factors_path);
    if (! in.good()) {
        fprintf(stderr, "File %s does not exist\n", factors_path);
        return 1;
    }

    FingerGraph* graph;
    print_time();
    fprintf(stderr, "Building graph...\n");
    graph = new FingerGraph(in, k, limit, !no_norm, !no_enriched);
    print_time();
    fprintf(stderr, "Done\n");
    print_time();
    fprintf(stderr, "Printing graph...\n");
    graph->save();
    print_time();
    fprintf(stderr, "Done\n");


    return 0;
}
