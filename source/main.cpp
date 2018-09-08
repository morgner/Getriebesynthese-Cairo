#include "canvas.h"
#include <gtkmm/application.h>
#include <gtkmm/window.h>


int main(int argc, char** argv)
{

    auto app = Gtk::Application::create(argc, argv, "org.gtkmm.3lsynth");

    Gtk::Window window;
    window.resize(800,600);
    window.set_title("3-Lagen Synthese");

    CCanvas area;
    window.add(area);
    area.show();

    return app->run(window);

    // ------------------------------------------------

}
