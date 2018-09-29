#include "canvas.h"
#include "grundpunkt.h"
#include <string>
#include <array>
#include <math.h>
#include <pangomm/fontdescription.h>


VEbenenLagen    g_vEbenenLagen;// E1,  E2,  E3  = 3 homologe Lagen einer Ebene
VPolDreieck     g_vPolDreieck; // P12, P13, P23 = 3 Polpunkte zu den o.g.Ebenenlagen
SUmkreis        g_tUmkreis;
VGelenke        g_vGelenke;
SCollision      g_tCollision{SCollision::EWhat::none, 0, 0};
CGrundpunkt     g_oGrundpunkt({-1, -1});
VGrundpunkte    g_vGrundPunkte;
double          g_dGVS;
bool            g_bShowText{true};
bool            g_bShowHints{true};
bool            g_bShowBlink{true};
double dScale  {1.0}; // scale
bool   bTransInitialized {false}; // translation initialized?
SPoint dTrans     {0.0, 0.0}; // translation relativ pos
SPoint dTransStart{0.0, 0.0}; // translation start pos
SPoint dBaseMove  {0.0, 0.0}; // base pos for moving the canvas


SPoint g_tPointA;
SPoint g_tPointB;
bool   g_bCSplit{false};

double dCollisonDistance{0.0}; // collision distance

SPoint dCollisionOffset{0.0, 0.0}; // distance to a specified point

SPoint dmos{0.0, 0.0}; // mouse offset scaled

auto CalcDistance(SPoint const & a, SPoint const & b)
    {
    return sqrt( pow((a.x-b.x),2) + pow((a.y-b.y),2) );
    }

auto CalcAlpha(double const & a, double const & b, double const & c)
    {
    auto cosinus = (pow(b,2)+pow(c,2)-pow(a,2)) / (2*b*c);
    if ((cosinus<-1)||(cosinus>1))
	{
//	std::cout << "cos(alfa): " << cosinus << ", ";
	return 0.0;
	}

//    if (s<0) return (2*M_PI) - (acos( cosinus ));
    return acos( cosinus );
    }

auto CalcVectorSlope(SPoint const & a, SPoint const & b)
    {
    auto const dba = SPoint{ b.x-a.x, b.y-a.y };
    auto       sba = atan(dba.y / ((abs(dba.x)<.00001)?.00001:dba.x));

    if      ( (dba.x>0) && (dba.y<0) )
	sba = -sba;
    else if ( (dba.x<0) && (dba.y<0) )
	sba = M_PI - sba;
    else if ( (dba.x<0) && (dba.y>0) )
	sba = M_PI - sba;
    else if ( (dba.x>0) && (dba.y>0) )
	sba = 2*M_PI - sba;

    return sba;
    }

auto CalcVectorDiff(SPoint const & a, SPoint const & b, SPoint const & c)
    {
    auto slb = CalcVectorSlope(a,b);
    auto slc = CalcVectorSlope(a,c);
    return slc - slb;
    }

double MouseDistance( SPoint const & p, SPoint const & pMouse )
    {
    dCollisionOffset.x = (p.x - pMouse.x)*dScale;
    dCollisionOffset.y = (p.y - pMouse.y)*dScale;
    return CalcDistance(p, pMouse);
    }

SCollision MouseCollision(SPoint const & crtMousePoint)
    {
    SCollision ac;
    ac.eWhat  = SCollision::EWhat::none;
    int cnt{0};

    dCollisonDistance = 12.0;

    if ( g_vGrundPunkte.size() > 1 )
	{
	SPoint const & A = g_vGrundPunkte[0].GPoint(0);
	if ( MouseDistance({g_tPointA.x, g_tPointA.y}, crtMousePoint) < dCollisonDistance )
	    {
	    ac.eWhat  = SCollision::EWhat::A;
	    ac.nIndex = cnt;
	    ac.nSubIx = 0;
	    dmos.x = dCollisionOffset.x;
	    dmos.y = dCollisionOffset.y;
	    return std::move(ac);
	    }
	}

    for (auto const & a:g_vEbenenLagen)
	{
	if ( MouseDistance({a.x1, a.y1}, crtMousePoint) < dCollisonDistance )
	    {
	    ac.eWhat  = SCollision::EWhat::Ebene;
	    ac.nIndex = cnt;
	    ac.nSubIx = 0;
	    dmos.x = dCollisionOffset.x;
	    dmos.y = dCollisionOffset.y;
	    break;
	    }
	if ( MouseDistance({a.x2, a.y2}, crtMousePoint) < dCollisonDistance )
	    {
	    ac.eWhat  = SCollision::EWhat::Ebene;
	    ac.nIndex = cnt;
	    ac.nSubIx = 1;
	    dmos.x = dCollisionOffset.x;
	    dmos.y = dCollisionOffset.y;
	    break;
	    }
	SPoint m = a.M();
	if ( MouseDistance({m.x, m.y}, crtMousePoint) < dCollisonDistance )
	    {
	    ac.eWhat  = SCollision::EWhat::Ebene;
	    ac.nIndex = cnt;
	    ac.nSubIx = 2;
	    dmos.x = dCollisionOffset.x;
	    dmos.y = dCollisionOffset.y;
	    break;
	    }
	++cnt;
	}
    cnt = 0;
    for (auto const & a:g_vGrundPunkte)
	{
	if ( MouseDistance( a.P123(), crtMousePoint) < dCollisonDistance )
	    {
	    ac.eWhat  = SCollision::EWhat::Grundpunkt;
	    ac.nIndex = cnt;
	    ac.nSubIx = 0;
	    dmos.x = dCollisionOffset.x;
	    dmos.y = dCollisionOffset.y;
	    break;
	    }
	++cnt;
	}

    return std::move(ac);
    }

void CCanvas::MovePunktA(SPoint const & tPointsTo)
    {
    if ( g_vGrundPunkte.size() < 1 ) return;

    SPoint const & A0 = g_vGrundPunkte[0].G0();
    auto w = CalcVectorSlope(A0, tPointsTo) + M_PI/2;
         m_tAnimator = w / (2*M_PI) ;
    }

void CCanvas::MoveEbenenPunkt(SPoint const & mp, double const & L)
    {
    SEbene & er = g_vEbenenLagen[g_tCollision.nIndex];

    SEbene tL{};
    if ( g_tCollision.nSubIx == 0 ) tL = {mp.x,mp.y,er.x2,er.y2};
    if ( g_tCollision.nSubIx == 1 ) tL = {er.x1,er.y1,mp.x,mp.y};
    if ( g_tCollision.nSubIx == 2 )
	{
	SPoint m = er.M();
	auto const dx = m.x - mp.x;
	auto const dy = m.y - mp.y;
	er = {er.x1-dx,er.y1-dy,er.x2-dx,er.y2-dy};
	}
    else
	{
	FixedLenLine(tL, L, g_tCollision.nSubIx == 0);
	er = {tL.x1,tL.y1,tL.x2,tL.y2};
	}
    }

void ExportSCAD( SPoint const & A0,
		 SPoint const & B0,
		 SPoint const & A,
		 SPoint const & B,
		 SEbene const & E1,
		 SEbene const & E2,
		 SEbene const & E3
		)
    {
    auto GL = CalcDistance(A0, B0);
    auto AL = CalcDistance(A0, A );
    auto BL = CalcDistance(B0, B );
    auto CL = CalcDistance(A,  B );

    TRenderItem oSubM{};
    TRenderData oData{};

    oSubM.emplace("A0x", std::to_string(A0.x));
    oSubM.emplace("A0y", std::to_string(A0.y));

    oSubM.emplace("B0x", std::to_string(B0.x));
    oSubM.emplace("B0y", std::to_string(B0.y));

    oSubM.emplace("GL",  std::to_string(GL));
    oSubM.emplace("AL",  std::to_string(AL));
    oSubM.emplace("BL",  std::to_string(BL));
    oSubM.emplace("CL",  std::to_string(CL));

    oData.emplace( "Getriebe", oSubM );
    oSubM.clear();

    oSubM.emplace("E1P1.x",  std::to_string(E1.x1));
    oSubM.emplace("E1P1.y",  std::to_string(E1.y1));
    oSubM.emplace("E1P2.x",  std::to_string(E1.x2));
    oSubM.emplace("E1P2.y",  std::to_string(E1.y2));

    oSubM.emplace("E2P1.x",  std::to_string(E2.x1));
    oSubM.emplace("E2P1.y",  std::to_string(E2.y1));
    oSubM.emplace("E2P2.x",  std::to_string(E2.x2));
    oSubM.emplace("E2P2.y",  std::to_string(E2.y2));

    oSubM.emplace("E3P1.x",  std::to_string(E3.x1));
    oSubM.emplace("E3P1.y",  std::to_string(E3.y1));
    oSubM.emplace("E3P2.x",  std::to_string(E3.x2));
    oSubM.emplace("E3P2.y",  std::to_string(E3.y2));

    oData.emplace( "Ebene", oSubM );
    oSubM.clear();

    Cte ote(oData, "3LagenSynthese.tmpl", "../templates/");
    std::ofstream f("./synthese.scad");
    f << ote << '\n';
    f.close();
//    std::cout << ote << '\n';
    }




bool CCanvas::on_key_press_event(GdkEventKey* key_event)
    {
    auto c = gdk_unicode_to_keyval (key_event->keyval);
    dScale = 1;
    dTrans = {0.0, 0.0};
    queue_draw();
    return true;
    }

bool CCanvas::on_button_press_event(GdkEventButton *event)
    {
    for ( auto const & a:m_voButtons )
	{
	if ( a.Collision(*event) )
	    {
	    m_oButtonPressed = a.text;
	    return true;
	    }
	}

    SPoint e    = (SPoint{*event} - dTrans)/dScale;
    dBaseMove   = *event;
    dTransStart = dTrans;

    if( (event->type == GDK_BUTTON_PRESS) && (event->button == 1) )
	{
        if ( !m_bFirstClick )
            {
            x1=x2=e.x;
            y1=y2=e.y;
            }
        if (g_tCollision.eWhat == SCollision::EWhat::none)
            {
	    switch (m_ePhase)
		{
		case EPhase::GrundPunkte:
		    g_oGrundpunkt.Update(e);
		    [[fallthrough]];
		case EPhase::EbenenLagen:
		    if ( !m_bFirstClick )
			{
			m_bFirstClick=true;
			}
		    break;
		default: break;
		}
            }
	}
    queue_draw();
    return true;
    }

SPoint g{0.0, 0.0};

bool CCanvas::on_motion_notify_event(GdkEventMotion *event)
    {
//    auto const e.x=event->x/dScale-dTrans.x/dScale;
    auto const e{ (SPoint{*event} - dTrans)/dScale };

    g = *event;
    if ( m_oButtonPressed.size() > 0 )
	{
	return true;
	}

    if ( event->type & GDK_MOTION_NOTIFY )
	if ( event->state & GDK_BUTTON1_MASK )
	    {
	    if ( m_bWithTraces && g_tCollision.eWhat != SCollision::EWhat::A )
	    	{
	    	m_vSpurE1.clear();
	    	m_vSpurE2.clear();
	    	m_tAnimStep = 0;
	    	}
	    switch (g_tCollision.eWhat)
		case SCollision::EWhat::none:
		    {
		    switch (m_ePhase)
			{
			case EPhase::EbenenLagen:
			    x2=e.x;
			    y2=e.y;
			    if ( m_nLenEbene != 0.0 )
				{
				SEbene tEbene{x1, y1, x2, y2};
				FixedLenLine(tEbene, m_nLenEbene);
				x2 = tEbene.x2;
				y2 = tEbene.y2;
				}
			    break;

			case EPhase::GrundPunkte:
			    g_oGrundpunkt.Update(e);
			    break;

			default:
			    dTrans = dTransStart - (dBaseMove - *event);
			    break;
			}
		    break;

		case SCollision::EWhat::A:
		    MovePunktA(e+dmos/dScale);
		    break;

		case SCollision::EWhat::Ebene:
		    MoveEbenenPunkt(e+dmos/dScale, m_nLenEbene);
		    break;

		case SCollision::EWhat::Grundpunkt:
		    g_vGrundPunkte[g_tCollision.nIndex].Update(e+dCollisionOffset/dScale);
		    break;
		}
	    }
	else
	    {
	    g_tCollision = MouseCollision(e);
	    }

    queue_draw(); // to be deleted
    return true;
    }

bool CCanvas::on_button_release_event(GdkEventButton* event)
    {
    if ( m_oButtonPressed.size() > 0 )
	{
	if ( m_oButtonPressed == "0" ) { m_dAnimate *= 0.9; m_dAnimate = (m_dAnimate<m_dAnimateMin)?m_dAnimateMin:m_dAnimate; }
	if ( m_oButtonPressed == "1" ) { m_dAnimate *= 1.1; m_dAnimate = (m_dAnimate>m_dAnimateMax)?m_dAnimateMax:m_dAnimate; }
	if ( m_oButtonPressed == "2" ) { m_bDurchschlagen = !m_bDurchschlagen; m_vSpurE1.clear(); m_vSpurE2.clear(); m_tAnimStep = 0;  }
	if ( m_oButtonPressed == "3" )   m_bDirectionLeft = !m_bDirectionLeft;
	if ( m_oButtonPressed == "4" )   m_bAnimate       = !m_bAnimate;
	if ( m_oButtonPressed == "5" ) { m_bWithTraces    = !m_bWithTraces; m_vSpurE1.clear(); m_vSpurE2.clear(); m_tAnimStep = 0; }
	if ( m_oButtonPressed == "6" )   m_bRotate        = !m_bRotate;
	if ( m_oButtonPressed == "7" )   g_bShowText      = !g_bShowText;
	if ( m_oButtonPressed == "8" )   g_bShowHints     = !g_bShowHints;
	if ( m_oButtonPressed == "9" )   g_bShowBlink     = !g_bShowBlink;
	if ( m_oButtonPressed =="11" )   m_bShowMouse     = !m_bShowMouse;

	m_oButtonPressed="";
	queue_draw();
	return true;
	}

    if ( (!m_bWithTraces) || (m_bWithTraces && g_tCollision.eWhat != SCollision::EWhat::A) )
	{
	m_vSpurE1.clear();
	m_vSpurE2.clear();
	m_tAnimStep = 0;
	}

    if (g_tCollision.eWhat == SCollision::EWhat::none)
	{
	if( (event->type == GDK_BUTTON_RELEASE) /*&& (event->button == 1)*/ )
	    {
	    switch (m_ePhase)
		{
		case EPhase::GrundPunkte:
		    g_oGrundpunkt.Update( {event->x/dScale-dTrans.x/dScale, event->y/dScale-dTrans.y/dScale} );
		    g_vGrundPunkte.emplace_back(g_oGrundpunkt);
		    if ( g_vGrundPunkte.size()  > 1 ) { m_ePhase = EPhase::Collision; }
		    break;

		case EPhase::EbenenLagen:
		    g_vEbenenLagen.emplace_back(SEbene{ x1, y1, x2, y2 });
		    x1 = y1 = x2 = y2 = 0;
		    if ( g_vEbenenLagen.size() == 1 )
			{
			double dx = g_vEbenenLagen[0].x1 - g_vEbenenLagen[0].x2;
			double dy = g_vEbenenLagen[0].y1 - g_vEbenenLagen[0].y2;
			m_nLenEbene = sqrt( dx*dx + dy*dy );
			}

		    switch ( g_vEbenenLagen.size() )
			{
			case 2: g_vPolDreieck.emplace_back( CalcPolpunkt(g_vEbenenLagen[1-1], g_vEbenenLagen[2-1]) );
				break;
			case 3: g_vPolDreieck.emplace_back( CalcPolpunkt(g_vEbenenLagen[1-1], g_vEbenenLagen[3-1]) );
				g_vPolDreieck.emplace_back( CalcPolpunkt(g_vEbenenLagen[2-1], g_vEbenenLagen[3-1]) );
				break;
			}
		    if ( g_vEbenenLagen.size()  > 2 ) { m_ePhase = EPhase::GrundPunkte; }
		    break;

		default: break;
		}
	    }
	}
    m_bFirstClick=false;

    queue_draw();

    if ( (g_vGrundPunkte.size()>1) && (g_vEbenenLagen.size()>2) )
	{
	SPoint A0 = g_vGrundPunkte[0].G0();
	SPoint B0 = g_vGrundPunkte[1].G0();
	SPoint A  = g_vGrundPunkte[0].GPoint(0);
	SPoint B  = g_vGrundPunkte[1].GPoint(0);
	ExportSCAD( A0, B0, A, B, g_vEbenenLagen[0], g_vEbenenLagen[1], g_vEbenenLagen[2] );
	}

    return true;
    }

bool CCanvas::on_scroll_event(GdkEventScroll *event)
    {
    Gtk::Allocation allocation = get_allocation();
    const int width  = allocation.get_width();
    const int height = allocation.get_height();

    SPoint const p0{ (*event - dTrans)/dScale };
    dScale *= (event->delta_y>0)?.9:1.1; if (dScale<.01) dScale=.01;
    SPoint const p1{ (*event - dTrans)/dScale };
    dTrans -= (p0-p1)*dScale;

    queue_draw();
    return true;
    }

template<typename P>
    void MoveTo(Cairo::RefPtr<Cairo::Context> const & cr, P const & tPoint)
	{
	cr->move_to(tPoint.x, tPoint.y);
	}

template<typename P>
    void LineTo(Cairo::RefPtr<Cairo::Context> const & cr, P const & tPoint)
	{
	cr->line_to(tPoint.x, tPoint.y);
	}

template<typename L>
    void Line(Cairo::RefPtr<Cairo::Context> const & cr, L const & tLine )
	{
	cr->move_to(tLine.x1,tLine.y1);
	cr->line_to(tLine.x2,tLine.y2);
	cr->stroke();
	}

template<typename S, typename... P>
    void Line(Cairo::RefPtr<Cairo::Context> const & cr, S const & tPoint1, P const &... tPoints )
	{
	 MoveTo(cr, tPoint1);
	(LineTo(cr, tPoints), ...);
	cr->stroke();
	}

template<typename S, typename... P>
    void Polygon(Cairo::RefPtr<Cairo::Context> const & cr, S const & tPoint1, P const &... tPoints )
	{
	 MoveTo(cr, tPoint1);
	(LineTo(cr, tPoints), ...);
	cr->fill();
	}

template<typename P>
    void Circle(Cairo::RefPtr<Cairo::Context> const & cr, P const & tPoint, double const & dRadius )
	{
	cr->arc(tPoint.x, tPoint.y, dRadius, 0, 2*M_PI);
	cr->fill();
	}

template<typename P>
    void Ring(Cairo::RefPtr<Cairo::Context> const & cr, P const & tPoint, double const & dRadius )
	{
	cr->arc(tPoint.x, tPoint.y, dRadius, 0, 2*M_PI);
	cr->stroke();
	}


void draw_text(Cairo::RefPtr<Cairo::Context> const & cr,
	       double posx, double posy,
	       std::string const & crsText, double const & scale = 1.0)
{
    if (!g_bShowText) return;

    cr->save();

//  PangoAttribute *attr = pango_attr_size_new_absolute(20 * PANGO_SCALE);

  // http://developer.gnome.org/pangomm/unstable/classPango_1_1FontDescription.html
  Pango::FontDescription font;

  font.set_family("Monospace");
  font.set_absolute_size(12.0 * PANGO_SCALE/scale);
//  font.set_weight(Pango::WEIGHT_BOLD);

  // http://developer.gnome.org/pangomm/unstable/classPango_1_1Layout.html
  CCanvas w;
  auto layout = w.create_pango_layout(crsText);

  layout->set_font_description(font);
  int text_width;
  int text_height;

  layout->get_pixel_size(text_width, text_height);

  cr->move_to(posx-text_width/2, posy-text_height/2);

  layout->show_in_cairo_context(cr);

  cr->restore();
}

void draw_text(Cairo::RefPtr<Cairo::Context> const & cr,
	       SPoint p,
	       std::string const & crsText, double const & scale = 1.0)
    {
    draw_text(cr, p.x, p.y, crsText, scale);
    }


void draw_ebene(Cairo::RefPtr<Cairo::Context> const & cr,
	        SEbene const & croEbene, int nId)
    {
    cr->set_line_cap(Cairo::LINE_CAP_ROUND);
    cr->set_line_width(13);
    cr->set_source_rgb(1,0,0);
    Line(cr, croEbene);

    cr->set_source_rgb(1,1,1);
    Circle(cr, SPoint{croEbene.x1,croEbene.y1},5);
    Circle(cr, SPoint{croEbene.x2,croEbene.y2},5);

    cr->set_source_rgb(0,0,0);
    draw_text(cr,  croEbene.x1+12, croEbene.y1-12, "p1");
    draw_text(cr,  croEbene.x2+12, croEbene.y2-12, "p2");
    draw_text(cr, (croEbene.x2+croEbene.x1)/2,
	          (croEbene.y2+croEbene.y1)/2, "E"+std::to_string(nId));
    }

void draw_grundpunkt(Cairo::RefPtr<Cairo::Context> const & cr,
	             CGrundpunkt & croGP, int nId)
    {
    croGP.Show( cr, g_vPolDreieck);

    SPoint const G0 = croGP.G0();

    std::string sId = std::array<std::string,5>{"U","A","B","C","E"}[(nId>4)?4:nId];

    if (g_bShowHints)
	{
	for ( int i{0}; i<3; ++i)
	    {
	    SPoint const GPoint = croGP.GPoint(i);

	    cr->set_line_cap(Cairo::LINE_CAP_ROUND);
	    cr->set_line_width(4);
	    cr->set_source_rgba(.5,.5,.5,.5);
	    Line(cr, G0, GPoint);

	    if (nId == 1)
		{
		if (g_vGrundPunkte.size()>1)
		    {
		    auto const A{g_vGrundPunkte[0].GPoint(i)};
		    auto const B{g_vGrundPunkte[1].GPoint(i)};
		    Line(cr, A, B);
		    if (i==0)
			{
			cr->set_source_rgb(0,0,0);
			draw_text(cr, (A+B)/2, std::to_string((int) CalcDistance(A,B)));
			}
		    }
		else
		    {
		    if (g_vGrundPunkte.size()>0)
			{
			auto const A{g_vGrundPunkte[0].GPoint(i)};
			auto const B{g_oGrundpunkt.GPoint(i)};
			Line(cr, A, B);
			if (i==0)
			    {
			    cr->set_source_rgb(0,0,0);
			    draw_text(cr, (A+B)/2, std::to_string((int) CalcDistance(A,B)));
			    }
			}
		    }
		}

	    if (i==0)
		{
		cr->set_source_rgb(0,0,0);
		draw_text(cr, (G0+GPoint)/2, std::to_string((int) CalcDistance(G0, GPoint) ));
		}
	    cr->set_source_rgba(.5,.5,.5,.5);
	    }
    //--------------------------------
	for ( int i{0}; i<3; ++i)
	    {
	    SPoint const GPoint = croGP.GPoint(i);

	    cr->set_source_rgb(.8,.8,1);
	    Circle(cr,GPoint,10);
	    cr->set_source_rgb(0,0,0);
	    cr->set_line_width(2);
	    Ring(cr,GPoint,10);
	    cr->set_source_rgb(0,0,0);
	    draw_text(cr, GPoint, sId+std::to_string(i+1));
	    }

	cr->set_source_rgb(0,0,0);
	draw_text(cr, croGP.P123().x,croGP.P123().y+27, sId+"123");
	}

    // G0 ------------------------------------------------
    cr->set_source_rgb(1,1,1);
    Circle(cr, G0,25);
    cr->set_source_rgb(0,0,0);
    cr->set_line_width(4);
    Ring(cr, G0,25);

    cr->set_line_width(1);
    cr->move_to(G0.x,G0.y);
    cr->line_to(G0.x+4,G0.y);
    cr->arc(G0.x,G0.y,25,0,M_PI/2);
    cr->line_to(G0.x,G0.y);
    cr->fill();
    cr->move_to(G0.x,G0.y);
    cr->line_to(G0.x-4,G0.y);
    cr->arc(G0.x,G0.y,25,M_PI,M_PI/2*3);
    cr->line_to(G0.x,G0.y);
    cr->fill();

    cr->set_source_rgb(0,0,0);
    draw_text(cr, G0.x,G0.y-32, sId+"0");
    }

void draw_pol(Cairo::RefPtr<Cairo::Context> const & cr,
	      SPoint const & pol)
    {
    if (!g_bShowHints) return;

    cr->set_line_width(1);
    cr->set_source_rgb(.9,.9,.9);
    Circle(cr, pol, 12);
    cr->set_source_rgb(0,0,0);
    Ring(cr, pol, 12);
    draw_text(cr, pol, "P12");
    }

void draw_poldreieck(Cairo::RefPtr<Cairo::Context> const & cr,
	             VPolDreieck const & croPD, int nId)
    {
    if (!g_bShowHints) return;

    cr->set_line_cap(Cairo::LINE_CAP_ROUND);
    cr->set_line_width(1);
    cr->set_source_rgba(0,1,0,.5);
    Polygon(cr, croPD[0], croPD[1], croPD[2]);
    cr->set_source_rgb(.5,.5,.5);
    Line(cr, croPD[0], croPD[1], croPD[2], croPD[0]);

    SUmkreis tUmkreis = Umkreis( croPD[0], croPD[1], croPD[2] );
    cr->set_source_rgba(1,1,1,.25);
    Circle(cr, tUmkreis.MidPnt, tUmkreis.Radius);
    cr->set_source_rgb(.5,.5,.5);
    Ring(cr, tUmkreis.MidPnt, tUmkreis.Radius);

    cr->set_source_rgb(0,0,0);
    draw_text(cr, croPD[0], "P12\n");
    draw_text(cr, croPD[1], "P13\n");
    draw_text(cr, croPD[2], "P23\n");
    }

bool CCanvas::on_draw(Cairo::RefPtr<Cairo::Context> const & cr)
    {
    Gtk::Allocation allocation = get_allocation();
    const int width  = allocation.get_width();
    const int height = allocation.get_height();

    static auto w0 = width/2;
    static auto h0 = height/2;

    if ( false == bTransInitialized )
	{
	w0 = dTrans.x = width/2;
	h0 = dTrans.y = height/2;
	bTransInitialized = true;
	}
    if ( (w0!=width/2) || (h0!=height/2) )
	{
	dTrans.x -= w0 - width/2;  w0 = width/2;
	dTrans.y -= h0 - height/2; h0 = height/2;
	}

    Cairo::Matrix matrix(1.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    matrix.scale(dScale,dScale);
    matrix.translate(dTrans.x/dScale, dTrans.y/dScale);

    if ( m_bRotate && g_vGrundPunkte.size()>1 )
	{
	SPoint const & A0 = g_vGrundPunkte[0].G0();
	SPoint const & B0 = g_vGrundPunkte[1].G0();
	           g_dGVS = CalcVectorSlope( A0, B0 );
	matrix.rotate(g_dGVS);
	}
    else
	{
	g_dGVS= 0;
	}

    cr->transform(matrix);

/*
    cr->set_source_rgb(.8, .8, 1);
    cr->rectangle(-150, -150, 300, 300);
    cr->fill();
*/

// static content
    if (g_vEbenenLagen.size() == 2)
	{
	draw_pol(cr, { CalcPolpunkt(g_vEbenenLagen[1-1], g_vEbenenLagen[2-1]) });
	}

    if (g_vEbenenLagen.size() == 3)
	{
	g_vPolDreieck[0] = ( CalcPolpunkt(g_vEbenenLagen[1-1], g_vEbenenLagen[2-1]) );
	g_vPolDreieck[1] = ( CalcPolpunkt(g_vEbenenLagen[1-1], g_vEbenenLagen[3-1]) );
	g_vPolDreieck[2] = ( CalcPolpunkt(g_vEbenenLagen[2-1], g_vEbenenLagen[3-1]) );
	draw_poldreieck(cr, g_vPolDreieck, 0);
	}

// interactive content
    if ( g_tCollision.eWhat == SCollision::EWhat::A )
	{
	cr->set_source_rgb(.7, .5, 0);
//	cr->set_line_width(11);
	cr->set_line_width(2);

	cr->arc(g_tPointA.x, g_tPointA.y, dCollisonDistance, 0, 2*M_PI);
	cr->fill();
	}

    if ( g_tCollision.eWhat == SCollision::EWhat::Ebene )
	{
	cr->set_source_rgb(.8,.6,1);
	cr->set_line_width(11);
	switch ( g_tCollision.nSubIx )
	    {
	    case 0:
		cr->arc(g_vEbenenLagen[g_tCollision.nIndex].x1,
			g_vEbenenLagen[g_tCollision.nIndex].y1,dCollisonDistance,0,2*M_PI);
		break;
	    case 1:
		cr->arc(g_vEbenenLagen[g_tCollision.nIndex].x2,
			g_vEbenenLagen[g_tCollision.nIndex].y2,dCollisonDistance,0,2*M_PI);
		break;
	    case 2:
		cr->set_line_cap(Cairo::LINE_CAP_ROUND);
		cr->set_line_width(dCollisonDistance*2);
		cr->move_to(g_vEbenenLagen[g_tCollision.nIndex].x1,g_vEbenenLagen[g_tCollision.nIndex].y1);
		cr->line_to(g_vEbenenLagen[g_tCollision.nIndex].x2,g_vEbenenLagen[g_tCollision.nIndex].y2);
		cr->stroke();
		break;
	    }
	cr->fill();
	}

    if ( g_tCollision.eWhat == SCollision::EWhat::Grundpunkt )
	{
	cr->set_source_rgb(.5,.5,0);
//	cr->set_line_width(11);
	cr->set_line_width(2);

	cr->arc(g_vGrundPunkte[g_tCollision.nIndex].P123().x,
		g_vGrundPunkte[g_tCollision.nIndex].P123().y, 21, 0, 2*M_PI);
	cr->fill();
	}

    int n=0;
    for (auto & a:g_vGrundPunkte)
	{
	draw_grundpunkt(cr, a, ++n);
	}

    n = 0;
    for (auto & a:g_vEbenenLagen)
	{
	draw_ebene(cr, a, ++n);
	}
    if (g_vEbenenLagen.size() == 2)
	{
	}


// dynamic content
    if ( m_bFirstClick )
	{
	switch (m_ePhase)
	    {
	    case EPhase::EbenenLagen:
		draw_ebene(cr, SEbene{x1,y1,x2,y2}, ++n);
		break;
	    case EPhase::GrundPunkte:
		if (m_bFirstClick)
		    {
		    draw_grundpunkt( cr, g_oGrundpunkt, g_vGrundPunkte.size()+1 );
		    }
		break;

	    default: break;
	    }
	}

/*
    cr->set_source_rgb(.8,.8,.0);
    cr->set_line_width(2);
    cr->arc(mx/s-tx/s, my/s-ty/s, 7/s,0,2*M_PI);
    cr->fill();

    cr->set_source_rgb(.0,.8,.8);
    cr->set_line_width(2);
    cr->arc(-tx/s, -ty/s, 11/s,0,2*M_PI);
    cr->fill();

    cr->set_source_rgb(0,0,0);
    cr->set_font_size(24.0/s);
    draw_text(cr, 200/s-tx/s, 50/s-ty/s, std::to_string((int)tx)   + ", " + std::to_string((int)ty));
    draw_text(cr, 200/s-tx/s, 65/s-ty/s, std::to_string((int)(tx/s)) + ", " + std::to_string((int)(ty/s)));
    draw_text(cr, 200/s-tx/s, 80/s-ty/s, std::to_string((int)mx)   + ", " + std::to_string((int)my));
    draw_text(cr, 200/s-tx/s, 95/s-ty/s, std::to_string((int)(mx/s)) + ", " + std::to_string((int)(my/s)));
*/
/*
    cr->set_source_rgb(.4, .4, 1);
    cr->rectangle(  -5,   -5, 10, 10);
    cr->fill();

    cr->save();
    cr->set_source_rgb(.8, .8, 1);
    cr->rotate(t*2*M_PI);
    cr->rectangle(  0,   0, 100, 10);
//    cr->rectangle(-50, -50, 100, 10);
    cr->fill();
    cr->restore();
*/
    if (g_vGrundPunkte.size()>1)
	{
	double lw=7;
	cr->set_line_width(lw);

	auto const S = sin(m_tAnimator*2*M_PI);
	auto const C = cos(m_tAnimator*2*M_PI);
	auto const & A0{g_vGrundPunkte[0].G0()};
	auto const & A1 {g_vGrundPunkte[0].GPoint(0)};
	auto const & A2 {g_vGrundPunkte[0].GPoint(1)};
	auto const & A3 {g_vGrundPunkte[0].GPoint(2)};
	auto const & B1 {g_vGrundPunkte[1].GPoint(0)};
	auto const & B2 {g_vGrundPunkte[1].GPoint(1)};
	auto const & B3 {g_vGrundPunkte[1].GPoint(2)};
	auto const & B0 {g_vGrundPunkte[1].G0()};

	auto const GL = CalcDistance( A0, B0 );
	auto const AL = CalcDistance( A0, A1 );
	auto const BL = CalcDistance( B0, B1 );
	auto const CL = CalcDistance( A1, B1 );
	auto const DL = CalcDistance( A1, B0 );

	    g_tPointA = SPoint{A0.x+AL*S, A0.y+AL*C};
	auto const  d = CalcDistance( g_tPointA , B0 );

	auto const epsi = CalcVectorSlope(g_tPointA, B0) -M_PI/2;

	auto const beta = CalcAlpha(BL,CL,d); // 143,  86, CalcDistance( g_tPointA , B0 )  g:58, a:43
	auto const gama = CalcAlpha(CL,BL,d); //  86, 143, CalcDistance( g_tPointA , B0 )

	std::cout << "e: " << epsi/M_PI*180 << ", b: " << beta/M_PI*180 << ", g: " << gama/M_PI*180 << '\n';

//	if ( beta*gama == 0.0 ) return true;

	SPoint g_tPointB;
	SPoint pab;
	if ( m_bDurchschlagen )
	    {
	    g_tPointB  = SPoint{B0.x+BL*sin(epsi-gama), B0.y+BL*cos(epsi-gama)};
	    pab = SPoint{g_tPointA.x+CL*sin(epsi+beta+M_PI), g_tPointA.y+CL*cos(epsi+beta+M_PI)};
	    }
	else
	    {
	    g_tPointB = SPoint{B0.x+BL*sin(epsi+gama), B0.y+BL*cos(epsi+gama)};
	    pab = SPoint{g_tPointA.x+CL*sin(epsi-beta+M_PI), g_tPointA.y+CL*cos(epsi-beta+M_PI)};
	    }
//	std::cout << epsi/M_PI*180 << ", a="  << (B0.y-pa.y) << ", b="  << (B0.x-pa.x) << ", c="  << d << ", " << '\n';

	cr->set_line_width(lw+2);
	cr->set_source_rgb(1,1,1);
	if ( d > (BL+CL) )
	    cr->set_source_rgb(0,0,0);

	if (g_bShowBlink)
	    {
	    Line(cr, A0, g_tPointA);
	    Line(cr, g_tPointA, pab);
	    Line(cr, B0, g_tPointB);
	    }
	else
	    {
	    cr->set_source_rgb(.5,.0,.0); Line(cr, A0, g_tPointA);
	    cr->set_source_rgb(.0,.5,.0); Line(cr, g_tPointA,pab);
	    cr->set_source_rgb(.0,.0,.5); Line(cr, B0, g_tPointB);
	    }

	if ( CalcDistance(A1, g_tPointA)<20 )
	    if ( CalcDistance(B1, g_tPointB)<30 )
		cr->set_source_rgb(.7, 1, .7);
	    else
		cr->set_source_rgb(1, .4, .7);
	else if ( CalcDistance(A2, g_tPointA)<20 )
	    if ( CalcDistance(B2, g_tPointB)<30 )
		cr->set_source_rgb(.7, 1, .7);
	    else
		cr->set_source_rgb(1, .4, .7);
	else if ( CalcDistance(A3, g_tPointA)<20 )
	    if ( CalcDistance(B3, g_tPointB)<30 )
		cr->set_source_rgb(.7, 1, .7);
	    else
		cr->set_source_rgb(1, .4, .7);
	else
	    cr->set_source_rgb(0,0,0);

	g_bCSplit = ( d > (BL+CL) );
	if ( g_bCSplit )
	    cr->set_source_rgb(1,0,0);

	cr->set_line_width(lw);

	if (g_bShowBlink)
	    {
	    Line(cr, A0, g_tPointA);
	    Line(cr, g_tPointA, pab);
	    Line(cr, B0, g_tPointB);
	    }
	else
	    {
	    cr->set_source_rgb(.5,.0,.0); Line(cr, A0, g_tPointA);
	    cr->set_source_rgb(.0,.5,.0); Line(cr, g_tPointA,pab);
	    cr->set_source_rgb(.0,.0,.5); Line(cr, B0, g_tPointB);
	    }

	Ring(cr, g_tPointA, 3);
	Ring(cr, g_tPointB, 3);

	cr->set_source_rgb(1,1,1);

	Circle(cr, g_tPointA, 3);
	Circle(cr, g_tPointB, 3);


	// Ebene **************************************
	SPoint const ep1{g_vEbenenLagen[0].x1,g_vEbenenLagen[0].y1};
	SPoint const ep2{g_vEbenenLagen[0].x2,g_vEbenenLagen[0].y2};
	auto       we  = CalcVectorSlope(ep1, ep2);
	auto       le  = CalcDistance(ep1, ep2);

	auto       wg  = CalcVectorSlope(g_tPointA, g_tPointB);
	auto       wb  = CalcVectorSlope(B1, ep1);
	auto       dw  = CalcVectorDiff(B1, A1, ep1);

	auto      dwe  = CalcVectorDiff(ep1, B1, ep2);

	auto const lve = CalcDistance(B1, ep1);

//	std::cout << " dw: " << dw/M_PI*180 << ",  wg: " << wg/M_PI*180 << ",  wb: " << wb/M_PI*180 << ",  lve: " << lve << '\n';

	SPoint be1{g_tPointB.x+lve*sin(wg+dw-M_PI/2), g_tPointB.y+lve*cos(wg+dw-M_PI/2)};
	SPoint be2{be1.x+le*sin(wg+dw+dwe+M_PI/2), be1.y+le*cos(wg+dw+dwe+M_PI/2)};
	cr->set_source_rgba(0,1,1,.4);
	Polygon(cr, g_tPointA, be1, g_tPointB);

	cr->set_source_rgba(.25,.75,.75,.75);
	cr->set_line_width(2);
	Line(cr, g_tPointA, be1, g_tPointB);

	cr->set_source_rgba(.99,.2,.2,.5);
	cr->set_line_width(13);
	Line(cr, be1, be2);


	if ( m_bWithTraces )
	    {
	    if (360*m_tAnimStep+1 > m_vSpurE1.size())
		{
		m_vSpurE1.emplace_back(SPointT{be1.x,be1.y,g_bCSplit});
		m_vSpurE2.emplace_back(SPointT{be2.x,be2.y,g_bCSplit});
		}
	    if (m_tAnimStep<=1.1) m_tAnimStep += m_dAnimate;
	    cr->set_source_rgb(0,0,0);

	    cr->set_line_width(2);
	    bool b;
	    if (m_vSpurE1.size())
		{
		cr->move_to(m_vSpurE1[0].x, m_vSpurE1[0].y);
		b = m_vSpurE1[0].bCSplit;
		}
	    for (auto const & a:m_vSpurE1)
		{
		if (b) cr->set_source_rgb(1,0,0); else cr->set_source_rgb(0,0,0);
		if ( a.bCSplit != b )
		    {
		    cr->line_to(a.x, a.y);
		    cr->stroke();
		    }
		cr->line_to(a.x, a.y);
		b = a.bCSplit;
		}
	    cr->stroke();

	    if (m_vSpurE2.size())
		{
		cr->move_to(m_vSpurE2[0].x, m_vSpurE2[0].y);
		b = m_vSpurE2[0].bCSplit;
		}
	    for (auto const & a:m_vSpurE2)
		{
		if (b) cr->set_source_rgb(1,0,0); else cr->set_source_rgb(0,0,0);
		if ( a.bCSplit != b )
		    {
		    cr->line_to(a.x, a.y);
		    cr->stroke();
		    }
		cr->line_to(a.x, a.y);
		b = a.bCSplit;
		}
	    cr->stroke();
	    }
	}



    cr->rotate(-g_dGVS);

    //  GUI
    auto constexpr uiOffset{ 32.0};
    auto constexpr uiBaseWd{  8.0};
    auto constexpr uiBaseLn{128.0};

    auto const dAniGui = (uiBaseLn-2*uiBaseWd)*m_dAnimate/(m_dAnimateMax-m_dAnimateMin);


    // Buttons
    int i{0}, nGrBtns{9};
    for ( auto const & a:m_voButtons )
	{
	++i;

	if ( a.Collision({g.x,g.y}) )
	    {
	    cr->set_source_rgb(0,1,0);
	    cr->rectangle( (a.x-dTrans.x)/dScale, (a.y-dTrans.y)/dScale, (a.w)/dScale, (a.h)/dScale );
	    cr->fill();
	    }


	static double lastScale{-1};
	static Glib::RefPtr<Gdk::Pixbuf> const image  = Gdk::Pixbuf::create_from_file("buttons.png");
	static Glib::RefPtr<Gdk::Pixbuf>       imageS = image->scale_simple( (a.w)/dScale*nGrBtns, (a.h)/dScale, Gdk::INTERP_BILINEAR);
	if ( lastScale != dScale )
	    {
	    imageS = image->scale_simple( (a.w)/dScale*nGrBtns, (a.h)/dScale, Gdk::INTERP_BILINEAR);
	    lastScale = dScale;
	    }
	if (i < nGrBtns+1)
	    {
	    Gdk::Cairo::set_source_pixbuf(cr, imageS, (a.x-a.w*(i-1)-dTrans.x)/dScale, (a.y-dTrans.y)/dScale);

	    cr->rectangle( (a.x-dTrans.x)/dScale, (a.y-dTrans.y)/dScale, (a.w)/dScale, (a.h)/dScale );
	    cr->fill();
	    }
	else
	    {
	    if ( !a.Collision({g.x,g.y}) )
		{
		cr->set_source_rgb(.8, .8, .9);
		cr->rectangle( (a.x-dTrans.x)/dScale, (a.y-dTrans.y)/dScale, (a.w)/dScale, (a.h)/dScale );
		cr->fill();
		}
	    cr->set_source_rgb(0,0,0);
	    draw_text(cr,  (a.x+a.w/2-dTrans.x)/dScale, (a.y+a.h/2-dTrans.y)/dScale, a.text, dScale);
	    }
	}

    cr->set_line_cap(Cairo::LINE_CAP_ROUND);
    cr->set_source_rgb(0,0,0);
    cr->set_line_width(uiBaseWd/dScale);

    cr->move_to((uiOffset-dTrans.x)/dScale, (uiOffset-dTrans.y         )/dScale);
    cr->line_to((uiOffset-dTrans.x)/dScale, (uiOffset-dTrans.y+uiBaseLn)/dScale);
    cr->stroke();


    cr->set_line_width(5.0/dScale);

    cr->move_to((uiOffset-dTrans.x           )/dScale, (uiOffset-dTrans.y+uiBaseWd/2)/dScale);
    cr->line_to((uiOffset-dTrans.x+2*uiBaseWd)/dScale, (uiOffset-dTrans.y+uiBaseWd/2)/dScale);
    cr->stroke();

    cr->move_to((uiOffset-dTrans.x           )/dScale, (uiOffset-dTrans.y-uiBaseWd/2+uiBaseLn)/dScale);
    cr->line_to((uiOffset-dTrans.x+2*uiBaseWd)/dScale, (uiOffset-dTrans.y-uiBaseWd/2+uiBaseLn)/dScale);
    cr->stroke();

    cr->set_source_rgb(0,1,0);
    cr->set_line_width(1.5/dScale);

    cr->move_to((uiOffset-dTrans.x           )/dScale, (uiOffset-dTrans.y-uiBaseWd/2+uiBaseLn-dAniGui)/dScale);
    cr->line_to((uiOffset-dTrans.x+2*uiBaseWd)/dScale, (uiOffset-dTrans.y-uiBaseWd/2+uiBaseLn-dAniGui)/dScale);
    cr->stroke();

    if ( m_bShowMouse )
	{
	cr->set_source_rgb(.75,.75,0);
	Circle(cr,g,3);
	}

    return true;
    }
