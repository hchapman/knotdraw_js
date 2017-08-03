#include <emscripten.h>

extern int EMSCRIPTEN_KEEPALIVE randomDiagram(int n_crossings,
                                              int n_components,
                                              int max_att,
                                              int dia_type,
                                              int seed,
                                              int32_t** vertData);
