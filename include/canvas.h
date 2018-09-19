#ifndef __CANVAS_H
#define __CANVAS_H

#include <gtkmm/drawingarea.h>

struct SCollision
    {
    enum class EWhat {
	none,
	Ebene,
	Grundpunkt,
    } eWhat;
    int nIndex; // 1,2,3 = E1, E2, E3 ; 1,2 =G1, G2
    int nSubIx; // 1,2,3 = P1, P2, Pm
    };

struct SPointB
    {
    SPointB() = default;
    template<typename T>
	SPointB(T x, T y) : x(x), y(y) {}
    double x{0}, y{0};
    };

struct SButton
    {
    double x,y,w,h;
    std::string t;

    SButton (double const & ix, double const & iy,
	     double const & iw, double const & ih,
	     std::string it)
	: x(ix), y(iy), w(iw), h(ih), t(it)
	{}
    bool Collision(SPointB const & p) const
	{
	return ( p.x > x && p.x < x+w && p.y > y && p.y < y+h );
	}
    };

using VButtons     = std::vector<SButton>;


class CCanvas : public Gtk::DrawingArea
{
    double           t{0},t0{0}; // $t animation parameter
    sigc::slot<bool> m_fSlot;
    sigc::connection m_fConnection;
    bool Animate(int);

    VButtons    m_voButtons;
    std::string m_oButtonPressed;
    bool        m_bDurchschlagen; // 1
    bool        m_bDirectionLeft;
    bool        m_bAnimate;

    std::vector<SPointB> m_vSpurE1;
    std::vector<SPointB> m_vSpurE2;

public:
  CCanvas();
  virtual ~CCanvas() {m_fConnection.disconnect();};
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
