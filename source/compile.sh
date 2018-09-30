gcc -I ../include/ `pkg-config --libs --cflags gtkmm-3.0 cairomm-1.0` -std=c++2a -lstdc++ -lm main.cpp canvas.cpp getrgeomath.cpp
