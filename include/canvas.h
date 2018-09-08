#ifndef __CANVAS_H
#define __CANVAS_H

#include <gtkmm/drawingarea.h>

struct SAllCollision
    {
    enum class EWhat {
	none,
	Ebene,
	Grundpunkt,
    } eWhat;
    int nIndex; // 1,2,3 = E1, E2, E3 ; 1,2 =G1, G2
    int nSubIx; // 1,2,3 = P1, P2, Pm
    };

class CCanvas : public Gtk::DrawingArea
{
    double  t{0}; // $t animation parameter

public:
  CCanvas();
  virtual ~CCanvas() {};
  void MoveEbenenPunkt(double const & x,double const & y,double const & L);

protected:
  //Override default signal handler:
  bool on_draw(const Cairo::RefPtr<Cairo::Context>& cr) override;
  bool on_scroll_event(GdkEventScroll *event) override;
  bool on_button_press_event(GdkEventButton * event) override;
  bool on_motion_notify_event(GdkEventMotion *event) override;
  bool on_button_release_event(GdkEventButton* release_event) override;
  bool on_key_press_event(GdkEventKey* key_event) override;

private:
  void draw_rectangle(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);


private:

  enum class EPhase
      {
      EbenenLagen,
      GrundPunkte,
      Collision,
      };
  EPhase m_ePhase = EPhase::EbenenLagen;

  double m_nLenEbene{0};

    //display Pixbuf
    Glib::RefPtr < Gdk::Pixbuf > display;

    //two coordinates
    double x1;
    double y1;
    double x2;
    double y2;

    bool m_bFirstClick;

};

#endif // __CANVAS_H
