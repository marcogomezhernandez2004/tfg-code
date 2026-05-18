PLUGIN_NAME = bidirectional_chemical_synapse_BO

HEADERS = bidirectional_chemical_synapse_BO.h\
          utils.hpp\

SOURCES = bidirectional_chemical_synapse_BO.cpp\
          moc_bidirectional_chemical_synapse_BO.cpp\
          BO.cpp\
          evaluation.cpp\

LIBS = -lkfr_dsp_avx2 -lboost_serialization -lboost_filesystem -lboost_system -lnlopt

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile

CXXFLAGS := $(filter-out -std=c++11,$(CXXFLAGS)) \
            -std=c++20 -O2 -march=native -DNDEBUG -DUSE_NLOPT \
            -Iextern/limbo/src $(shell pkg-config --cflags eigen3) \
            -I$(HOME)/.local/include -Iinclude
