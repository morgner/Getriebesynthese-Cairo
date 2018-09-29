#ifndef __CANVAS_H
#define __CANVAS_H

#include "getrgeomath.h"

// #include <string>
#include <gtkmm.h>
#include <gtkmm/drawingarea.h>




struct SPointB
    {
    SPointB() = default;
    template<typename T>
	SPointB(T x, T y, bool b=false) : x(x), y(y), bCSplit(b)  {}
    template<typename T>
	SPointB(T const & t) : x(t.x), y(t.y) {}
    double x{0}, y{0};
    bool const bCSplit{false};
    };

struct SButton
    {
    double x,y,w,h;
    std::string text;

    SButton (double const & ix, double const & iy,
	     double const & iw, double const & ih,
	     std::string it)
	: x(ix), y(iy), w(iw), h(ih), text(it)
	{}
    bool Collision(SPointB const & p) const
	{
	return ( p.x > x && p.x < x+w && p.y > y && p.y < y+h );
	}
    };

using VButtons     = std::vector<SButton>;


class CCanvas : public Gtk::DrawingArea
    {
    public:
    CCanvas()
        {
        add_events(Gdk::BUTTON_PRESS_MASK | Gdk::SCROLL_MASK | Gdk::SMOOTH_SCROLL_MASK);
        add_events(Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK);
        add_events(Gdk::BUTTON1_MOTION_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::POINTER_MOTION_MASK);
        add_events(Gdk::KEY_PRESS_MASK | Gdk::KEY_RELEASE_MASK);

        m_fSlot       = sigc::bind(sigc::mem_fun(*this, &CCanvas::Animate), 0);
        m_fConnection = Glib::signal_timeout().connect(m_fSlot, 40);

        for ( int i{0}; i<11; ++i)
	    {
	    auto constexpr bs{38.0};
	    auto constexpr uix{20.0},uiy{20.0},uiw{bs},uih{bs};
	    auto constexpr bo{8.0};
	    m_voButtons.emplace_back( 72+uix+i*(uiw+bo), uiy, uiw, uih, std::to_string(i) );
	    }
        }


	virtual ~CCanvas() {m_fConnection.disconnect();};
	void MoveEbenenPunkt(double const & x,double const & y,double const & L);
	void MovePunktA(double const & x,double const & y);

    protected:
	//Override default signal handler:
	bool on_draw(const Cairo::RefPtr<Cairo::Context>& cr) override;
	bool on_scroll_event(GdkEventScroll *event) override;
	bool on_button_press_event(GdkEventButton * event) override;
	bool on_motion_notify_event(GdkEventMotion *event) override;
	bool on_button_release_event(GdkEventButton* release_event) override;
	bool on_key_press_event(GdkEventKey* key_event) override;

    private:
	double      m_nLenEbene{0};

	bool        m_bFirstClick{false};

	VButtons    m_voButtons;	// button bar
	std::string m_oButtonPressed;
	bool        m_bDurchschlagen{false};
	bool        m_bDirectionLeft{true};
	bool        m_bAnimate{true};
	bool        m_bWithTraces{true};
	bool        m_bRotate{false};
	bool        m_bShowMouse{false};

	// animation clock
	double           t{0},t0{0}; // $t animation parameter
	sigc::slot<bool> m_fSlot;
	sigc::connection m_fConnection;

	double       m_dAnimate   {0.0025};
	double const m_dAnimateMax{0.025 };
	double const m_dAnimateMin{0.0025};
	bool Animate(int c)
	    {
	    if (!m_bAnimate) return true;
	    if (m_bDirectionLeft)
		t = (t <=  m_dAnimate) ? 1 : t-m_dAnimate;
	    else
		t = (t >=1-m_dAnimate) ? 0 : t+m_dAnimate;
	    queue_draw();
	    return true;
	    }
	// koppelkurveb
	std::vector<SPointB> m_vSpurE1;
	std::vector<SPointB> m_vSpurE2;

	/// mouse parameters
	enum class EPhase
	    {
	    EbenenLagen,
	    GrundPunkte,
	    Collision,
	    } m_ePhase {EPhase::EbenenLagen};

	double x1{0};
	double y1{0};
	double x2{0};
	double y2{0};


    }; // class CCanvas

// __CANVAS_H
#endif
