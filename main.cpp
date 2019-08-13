//
//  main.cpp
//  dna
//
//  Created by Jacob Mathai on 7/15/19.
//  Copyright © 2019 Jacob Mathai. All rights reserved.
//

#include "strand.h"

int main(int argc, const char * argv[])
{
    map<string, string> lookup = create_lookup();
    strand s(atoi(argv[1]), "DNA", lookup);
    cout << s << endl << endl;
    s.print_all_data(lookup);
}
