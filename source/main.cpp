#include "canvas.h"
#include <gtkmm/application.h>
#include <gtkmm/window.h>

static gboolean time_handler(GtkWidget *widget)
{
  gtk_widget_queue_draw(widget);

  return TRUE;
}

int main(int argc, char** argv)
{
   auto app = Gtk::Application::create(argc, argv, "org.gtkmm.example");

   Gtk::Window win;
   win.resize(800,600);
   win.set_title("3-Lagen Synthese");

   CCanvas area;
   win.add(area);
   area.show();

   return app->run(win);
}

