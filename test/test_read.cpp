#include "../Grid.h"

int main() {
    HistogramValue h;
    h.readFromFile("a.dat");
    h.writeToFile("b.dat");
    return 0;
}
