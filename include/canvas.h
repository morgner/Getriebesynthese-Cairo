#ifndef __CANVAS_H
#define __CANVAS_H

#include <gtkmm/drawingarea.h>

struct SCollision
    {
    enum class EWhat {
	none,
	A,
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
	SPointB(T x, T y, bool b=false) : x(x), y(y), bCSplit(b)  {}
    double x{0}, y{0};
    bool   bCSplit{false};
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
    public:
	CCanvas();
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

	bool        m_bFirstClick{false};

	VButtons    m_voButtons;	// button bar
	std::string m_oButtonPressed;
	bool        m_bDurchschlagen{false};
	bool        m_bDirectionLeft{true};
	bool        m_bAnimate{true};
	bool        m_bWithTraces{true};
	bool        m_bRotate{false};

	// animation clock
	double           t{0},t0{0}; // $t animation parameter
	sigc::slot<bool> m_fSlot;
	sigc::connection m_fConnection;

	double       g_dAnimate   {0.0025};
	double const g_dAnimateMax{0.025 };
	double const g_dAnimateMin{0.0025};
	bool Animate(int c)
	    {
	    if (!m_bAnimate) return true;
	    if (m_bDirectionLeft)
		t = (t<=  g_dAnimate) ? 1 : t-g_dAnimate;
	    else
		t = (t>=1-g_dAnimate) ? 0 : t+g_dAnimate;
	    queue_draw();
	    return true;
	    }

	std::vector<SPointB> m_vSpurE1;
	std::vector<SPointB> m_vSpurE2;

	enum class EPhase
	    {
	    EbenenLagen,
	    GrundPunkte,
	    Collision,
	    } m_ePhase = EPhase::EbenenLagen;

	double m_nLenEbene{0};

	//two coordinates
	double x1{0};
	double y1{0};
	double x2{0};
	double y2{0};


    }; // class CCanvas

// __CANVAS_H
#endif
