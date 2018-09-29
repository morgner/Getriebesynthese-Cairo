#ifndef __GRUNDPUNKT_H
#define __GRUNDPUNKT_H

#include "getrgeomath.h"

#include <vector>
#include <array>
#include <math.h>
#include "te.h"

#include <gtkmm/drawingarea.h>

using namespace std::string_literals;

using TRenderer = Cairo::RefPtr<Cairo::Context>;


/*
 * CGrundpunkt
 *
 *
 */
class CGrundpunkt
    {
    protected:

	double    m_dX{};
	double    m_dY{};
	SUmkreis  m_tUK{};
	A3Gelenke m_a3Gelenke{}; // Gelenkpunkte
	bool      m_bFixed{true};


    public:

	CGrundpunkt(SPoint const & P123) : m_dX(P123.x), m_dY(P123.y) {}
	CGrundpunkt(CGrundpunkt const & src) = default;

	void FixIt( bool const bFixit = true ) { m_bFixed = bFixit; }
	bool Isfix() const { return m_bFixed; }

	SPoint const   P123() const { return {m_dX, m_dY}; }
	SPoint const & G0() const { return m_tUK.MidPnt; }
	SPoint const & GPoint( int const & i ) const { return m_a3Gelenke[i]; }


	void UpdateAndShow( TRenderer const & cr, SPoint const & P123, VPolDreieck const & Poldreieck)
	    {
	    Update(P123);
	    Show(cr, Poldreieck);
	    }

	void Update(SPoint const & P123)
	    {
	    m_dX = P123.x;
	    m_dY = P123.y;
	    }

	void Show( TRenderer const & cr, VPolDreieck const & Poldreieck)
	    {
	    cr->set_line_width(1); //!!

	    // P123
	    cr->set_source_rgba(0, .5, .5, 0.75);
	    cr->arc(m_dX, m_dY, 8, 0, 2*M_PI);
	    cr->fill();
	    cr->set_source_rgb (0, 0, 0);
	    cr->arc(m_dX, m_dY, 8, 0, 2*M_PI);
	    cr->stroke();

	    A3Gelenke & G = m_a3Gelenke;
	    for ( int n=0, i=0, j=0; n < 3; ++n )
		{
		// Lines: P12-P13, P12-P23, P13-P23
		switch (n)
		    {
		    case 0: i=0; j=1; break;
		    case 1: i=0; j=2; break;
		    case 2: i=1; j=2; break;
		    }
		SLine const PL{ (double)Poldreieck[i].x, (double)Poldreieck[i].y,
			        (double)Poldreieck[j].x, (double)Poldreieck[j].y};
		G[n] = PointMirror( { m_dX, m_dY }, PL );
		cr->set_source_rgba(1, 0, 0, 0.25);
		cr->arc(m_dX, m_dY, 20, 0, 2*M_PI);
		cr->fill();
		cr->arc(m_dX, m_dY, 20, 0, 2*M_PI);
		cr->set_source_rgb(0, 0, 0);
		cr->stroke();
		} // for (...

	    m_tUK = Umkreis( G[0], G[1], G[2] );
/*
	    cr->set_source_rgba(.5, .5, 0, 0.5);
	    cr->arc(m_tUK.MidPnt.x, m_tUK.MidPnt.y, 36, 0, 2*M_PI);
	    cr->fill();
	    cr->set_line_width(1);
	    cr->set_source_rgb(0, 0, 0);
	    cr->arc(m_tUK.MidPnt.x, m_tUK.MidPnt.y, 36, 0, 2*M_PI);
	    cr->stroke();
*/
	    }
    }; // class CGrundpunkt

using VGrundpunkte  = std::vector<CGrundpunkt>;
using A2Grundpunkte = std::array<CGrundpunkt, 2>;

// __GRUNDPUNKT_H
#endif
