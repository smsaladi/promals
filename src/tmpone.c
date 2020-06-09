#include "header_cpp.h"

class onealn {

        public:
                onealn();
                ~onealn();
                onealn *thisone;
};

onealn::onealn() {
        thisone = this;
}

onealn::~onealn() {
        ;
}

