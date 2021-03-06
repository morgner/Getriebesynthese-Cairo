#ifndef __CANVAS_H
#define __CANVAS_H

#include "getrgeomath.h"

#include <gtkmm.h>
#include <gtkmm/drawingarea.h>


struct SButton
    {
    double x,y,w,h;
    std::string const text;

    SButton (double const & ix, double const & iy,
	     double const & iw, double const & ih,
	     std::string const  s)
	: x(ix), y(iy), w(iw), h(ih), text(s)
	{}

    bool Collision(SPoint const & p) const
	{
	return ( p.x > x && p.x < x+w && p.y > y && p.y < y+h );
	}
    };

using VButtons = std::vector<SButton>;


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

	virtual ~CCanvas() { m_fConnection.disconnect(); };

	void MoveEbenenPunkt(SPoint const & mp, double const & L);
	void MovePunktA(SPoint const & tPointsTo);

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
	std::string m_oButtonPressed{""};
	bool        m_bDurchschlagen{false};
	bool        m_bDirectionLeft{true};
	bool        m_bAnimate{true};
	bool        m_bWithTraces{true};
	bool        m_bRotate{false};
	bool        m_bShowMouse{false};

	// animation clock
	double           m_tAnimator{0}; // $m_tAnimator animation parameter
	double           m_tAnimStep{0}; // intermediate animation parameter
	sigc::slot<bool> m_fSlot;
	sigc::connection m_fConnection;

	double           m_dAnimate   {0.0025}; // animation steps width
	double const     m_dAnimateMax{0.0250}; // minimal animation step width
	double const     m_dAnimateMin{0.0025}; // maximal animation step width
	bool Animate(int c)
	    {
	    if (!m_bAnimate) return true;
	    if (m_bDirectionLeft)
		m_tAnimator = (m_tAnimator <=  m_dAnimate) ? 1 : m_tAnimator-m_dAnimate;
	    else
		m_tAnimator = (m_tAnimator >=1-m_dAnimate) ? 0 : m_tAnimator+m_dAnimate;
	    queue_draw();
	    return true;
	    }

	// koppelkurven
	std::vector<SPointT> m_vSpurE1;
	std::vector<SPointT> m_vSpurE2;

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
