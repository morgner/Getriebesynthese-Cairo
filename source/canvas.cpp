#include "canvas.h"
#include "grundpunkt.h"
#include <string>
#include <array>
#include <math.h>
#include <gtkmm.h>
#include <pangomm/fontdescription.h>


VEbenenLagen    g_vEbenenLagen;// E1,  E2,  E3  = 3 homologe Lagen einer Ebene
VPolDreieck     g_vPolDreieck; // P12, P13, P23 = 3 Polpunkte zu den o.g.Ebenenlagen
SUmkreis        g_tUmkreis;
VGelenke        g_vGelenke;
SCollision      g_tCollision{SCollision::EWhat::none, 0, 0};
CGrundpunkt     g_oGrundpunkt({-1, -1});
VGrundpunkte    g_vGrundPunkte;

double dScale  {1.0}; // scale
bool   bTransInitialized {false}; // translation initialized?
double dTransX {0.0}; // translation relativ pos
double dTransY {0.0}; //
double dTransXStart{0.0}; // translation start pos
double dTransYStart{0.0}; //
double dBaseMoveX {0.0}; // base pos for moving the canvas
double dBaseMoveY {0.0}; //

double dCollisonDistance{0.0}; // collision distance

double dCollisionOffsetX{0}; // distance to a specified point
double dCollisionOffsetY{0};

double dmosx{0}; // mouse offset scaled
double dmosy{0};

template<typename A, typename B>
    auto CalcDistance(A const & a, B const & b)
	{
	return sqrt( pow((a.x-b.x),2) + pow((a.y-b.y),2) );
	}

template<typename A, typename B>
    auto CalcDistanceAndSign(A const & a, B const & b)
	{
	auto dbax = (b.x - a.x);
	auto dbay = (b.y - a.y);
	auto dbad = ( 0 > dbax*dbay ) ? -1 : 1;
	return dbad * sqrt( pow((a.x-b.x),2) + pow((a.y-b.y),2) );
	}


auto CalcAlpha(double const & a, double const & b, double const & c)
    {
    auto cosinus = (pow(b,2)+pow(c,2)-pow(a,2)) / (2*b*c);
    if ((cosinus<-1)||(cosinus>1))
	return 0.0;

    if (a<0) return (2*M_PI) - (acos( cosinus ));
    return acos( cosinus );
    }

auto CalcAlphaFull(double const & a, double const & b, double const & c)
    {
    auto cosinus = (pow(b,2)+pow(c,2)-pow(a,2)) / (2*b*c);
    if ((cosinus<-1)||(cosinus>1))
	return 0.0;

    if (a<0) return (2*M_PI) - (acos( cosinus ));
    return acos( cosinus );
    }

auto CalcVectorSlope(SPoint const & a, SPoint const & b)
    {
    auto const dba = SPoint{ b.x-a.x, b.y-a.y };
    auto       sba = atan(dba.y / ((abs(dba.x)<.0001)?.0001:dba.x));
//    if (dba.y<0) sba = 2*M_PI+sba;
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
/*
    auto const dba = SPoint{ b.x-a.x, b.y-a.y };
    auto       sba = atan(dba.y / ((abs(dba.x)<.0001)?.0001:dba.x));
    if (dba.y<0) sba = 2*M_PI+sba;

    auto const dca = SPoint{ c.x-a.x, c.y-a.y };
    auto       sca = atan(dca.y / ((abs(dca.x)<.0001)?.0001:dca.x));
    if (dba.y<0) sca = 2*M_PI+sca;

    return sca-sba;
*/
    auto slb = CalcVectorSlope(a,b);
    auto slc = CalcVectorSlope(a,c);
//    std::cout << ", a to ep1: " << slc/M_PI*180 << ", a to b: " << slb/M_PI*180 << " - ";
/*
    while (slb < 0)
	{
	slb += 2*M_PI;
	slc += 2*M_PI;
	}
*/
    return slc - slb;
    }

double MouseDistance( SPoint const & p, SPoint const & pMouse )
    {
    dCollisionOffsetX = (p.x - pMouse.x)*dScale;
    dCollisionOffsetY = (p.y - pMouse.y)*dScale;
    return CalcDistance(p, pMouse);
    }

SCollision MouseCollision(SPoint const & crtMousePoint)
    {
    SCollision ac;
    ac.eWhat  = SCollision::EWhat::none;
    int cnt{0};

    dCollisonDistance = 12.0;
    for (auto const & a:g_vEbenenLagen)
	{
	if ( MouseDistance({a.x1, a.y1}, crtMousePoint) < dCollisonDistance )
	    {
	    ac.eWhat  = SCollision::EWhat::Ebene;
	    ac.nIndex = cnt;
	    ac.nSubIx = 0;
	    dmosx = dCollisionOffsetX;
	    dmosy = dCollisionOffsetY;
	    break;
	    }
	if ( MouseDistance({a.x2, a.y2}, crtMousePoint) < dCollisonDistance )
	    {
	    ac.eWhat  = SCollision::EWhat::Ebene;
	    ac.nIndex = cnt;
	    ac.nSubIx = 1;
	    dmosx = dCollisionOffsetX;
	    dmosy = dCollisionOffsetY;
	    break;
	    }
	SPoint m = a.M();
	if ( MouseDistance({m.x, m.y}, crtMousePoint) < dCollisonDistance )
	    {
	    ac.eWhat  = SCollision::EWhat::Ebene;
	    ac.nIndex = cnt;
	    ac.nSubIx = 2;
	    dmosx = dCollisionOffsetX;
	    dmosy = dCollisionOffsetY;
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
	    dmosx = dCollisionOffsetX;
	    dmosy = dCollisionOffsetY;
	    break;
	    }
	++cnt;
	}

    return std::move(ac);
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

double       g_dAnimate   {0.0025};
double const g_dAnimateMax{0.025 };
double const g_dAnimateMin{0.0025};
bool CCanvas::Animate(int c)
    {
    if (!m_bAnimate) return true;
    if (m_bDirectionLeft)
	t = (t<=1-g_dAnimate) ? t-g_dAnimate : 1;
    else
	t = (t<=1-g_dAnimate) ? t+g_dAnimate : 0;
    queue_draw();
    return true;
    }

CCanvas::CCanvas()
    : m_bFirstClick(false),
      x1(0),y1(0),
      x2(0),y2(0),
      m_bDurchschlagen(false),
      m_bDirectionLeft(true),
      m_bAnimate(true)
    {
    add_events(Gdk::BUTTON_PRESS_MASK | Gdk::SCROLL_MASK | Gdk::SMOOTH_SCROLL_MASK);
    add_events(Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK);
    add_events(Gdk::BUTTON1_MOTION_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::POINTER_MOTION_MASK);
    add_events(Gdk::KEY_PRESS_MASK | Gdk::KEY_RELEASE_MASK);

    m_fSlot       = sigc::bind(sigc::mem_fun(*this, &CCanvas::Animate), 0);
    m_fConnection = Glib::signal_timeout().connect(m_fSlot, 40);


    for ( int i{0}; i<16; ++i)
	{
	auto constexpr uix{20.0},uiy{20.0},uiw{20.0},uih{20.0};
	auto constexpr bo{4.0};
	m_voButtons.emplace_back( 72+i*(uix+bo), uiy, uiw, uih, std::to_string(i) );
	}

    }


bool CCanvas::on_key_press_event(GdkEventKey* key_event)
    {
    auto c = gdk_unicode_to_keyval (key_event->keyval);
    dScale = 1;
    dTransX = dTransY = 0;
    queue_draw();
    return true;
    }

bool CCanvas::on_button_press_event(GdkEventButton *event)
    {
    for ( auto const & a:m_voButtons )
	{
	if ( a.Collision({event->x,event->y}) )
	    {
	    m_oButtonPressed=a.t;
	    return true;
	    }
	}

    auto const ex = event->x/dScale-dTransX/dScale;
               dBaseMoveX = event->x;
               dTransXStart = dTransX;
    auto const ey = event->y/dScale-dTransY/dScale;
               dBaseMoveY = event->y;
               dTransYStart = dTransY;

    if( (event->type == GDK_BUTTON_PRESS) && (event->button == 1) )
	{
        if ( !m_bFirstClick )
            {
            x1=x2=ex;
            y1=y2=ey;
            }
        if (g_tCollision.eWhat == SCollision::EWhat::none)
            {
	    switch (m_ePhase)
		{
		case EPhase::GrundPunkte:
		    g_oGrundpunkt.Update({ex, ey});
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

double gx{0};
double gy{0};
bool CCanvas::on_motion_notify_event(GdkEventMotion *event)
    {
    auto const ex=event->x/dScale-dTransX/dScale;
    auto const ey=event->y/dScale-dTransY/dScale;

    gx=event->x;
    gy=event->y;
    if ( m_oButtonPressed.size() > 0 )
	{
	return true;
	}

    if ( event->type & GDK_MOTION_NOTIFY )
	if ( event->state & GDK_BUTTON1_MASK )
	    {
	    switch (g_tCollision.eWhat)
		case SCollision::EWhat::none:
		    {
		    switch (m_ePhase)
			{
			case EPhase::EbenenLagen:
			    x2=ex;
			    y2=ey;
			    if ( m_nLenEbene != 0.0 )
				{
				SEbene tEbene{x1, y1, x2, y2};
				FixedLenLine(tEbene, m_nLenEbene);
				x2 = tEbene.x2;
				y2 = tEbene.y2;
				}
			    break;

			case EPhase::GrundPunkte:
			    g_oGrundpunkt.Update({ex, ey});
			    break;

			default:
			    dTransX = dTransXStart - (dBaseMoveX-event->x);
			    dTransY = dTransYStart - (dBaseMoveY-event->y);
			    break;
			}
		    break;

		case SCollision::EWhat::Ebene:
		    MoveEbenenPunkt(ex+dmosx/dScale, ey+dmosy/dScale, m_nLenEbene);
		    break;

		case SCollision::EWhat::Grundpunkt:
		    g_vGrundPunkte[g_tCollision.nIndex].Update({ex+dCollisionOffsetX/dScale, ey+dCollisionOffsetY/dScale});
		    break;
		}
	    }
	else
	    {
	    g_tCollision = MouseCollision({ex, ey});
	    }

    queue_draw(); // to be deleted
    return true;
    }

bool CCanvas::on_button_release_event(GdkEventButton* event)
    {
    if ( m_oButtonPressed.size() > 0 )
	{
//	std::cout << m_oButtonPressed << '\n';

	if ( m_oButtonPressed == "0" ) { g_dAnimate *= 0.9; g_dAnimate = (g_dAnimate<g_dAnimateMin)?g_dAnimateMin:g_dAnimate; }
	if ( m_oButtonPressed == "1" ) { g_dAnimate *= 1.1; g_dAnimate = (g_dAnimate>g_dAnimateMax)?g_dAnimateMax:g_dAnimate; }
	if ( m_oButtonPressed == "2" ) m_bDurchschlagen = !m_bDurchschlagen;
	if ( m_oButtonPressed == "3" ) m_bDirectionLeft = !m_bDirectionLeft;
	if ( m_oButtonPressed == "4" ) m_bAnimate       = !m_bAnimate;
	m_oButtonPressed="";
	return true;
	}

    if (g_tCollision.eWhat == SCollision::EWhat::none)
	{
	if( (event->type == GDK_BUTTON_RELEASE) /*&& (event->button == 1)*/ )
	    {
	    switch (m_ePhase)
		{
		case EPhase::GrundPunkte:
		    g_oGrundpunkt.Update( {event->x/dScale-dTransX/dScale, event->y/dScale-dTransY/dScale} );
		    g_vGrundPunkte.emplace_back(g_oGrundpunkt);
		    if ( g_vGrundPunkte.size()  > 1 ) { m_ePhase = EPhase::Collision; }
		    break;

		case EPhase::EbenenLagen:
		    g_vEbenenLagen.emplace_back(SEbene{ x1, y1, x2, y2 });
		    x1 = y1 = x2 = y2 = 0;
		    if ( g_vEbenenLagen.size() == 1 )
			{
			int dx = g_vEbenenLagen[0].x1 - g_vEbenenLagen[0].x2;
			int dy = g_vEbenenLagen[0].y1 - g_vEbenenLagen[0].y2;
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

    SPoint p0{event->x/dScale-dTransX/dScale, event->y/dScale-dTransY/dScale};
    dScale *= (event->delta_y>0)?.9:1.1; if (dScale<.01) dScale=.01;
    SPoint p1{event->x/dScale-dTransX/dScale, event->y/dScale-dTransY/dScale};
    dTransX -= (p0.x-p1.x)*dScale;
    dTransY -= (p0.y-p1.y)*dScale;

    queue_draw();
    return true;
    }

void CCanvas::MoveEbenenPunkt(double const & x,double const & y,double const & L)
    {
    SEbene & er = g_vEbenenLagen[g_tCollision.nIndex];

    SEbene tL{};
    if ( g_tCollision.nSubIx == 0 ) tL = {x,y,er.x2,er.y2};
    if ( g_tCollision.nSubIx == 1 ) tL = {er.x1,er.y1,x,y};
    if ( g_tCollision.nSubIx == 2 )
	{
	SPoint m = er.M();
	auto dx = m.x - x;
	auto dy = m.y - y;
	er = {er.x1-dx,er.y1-dy,er.x2-dx,er.y2-dy};
	}
    else
	{
	FixedLenLine(tL, L, g_tCollision.nSubIx == 0);
	er = {tL.x1,tL.y1,tL.x2,tL.y2};
	}
    }

void draw_text(Cairo::RefPtr<Cairo::Context> const & cr,
	       int posx, int posy,
	       std::string const & crsText, double const & scale = 1.0)
{
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

void draw_ebene(Cairo::RefPtr<Cairo::Context> const & cr,
	        SEbene const & croEbene, int nId)
    {
    cr->set_line_cap(Cairo::LINE_CAP_ROUND);
    cr->set_line_width(13);
    cr->set_source_rgb(1,0,0);
    cr->move_to(croEbene.x1,croEbene.y1);
    cr->line_to(croEbene.x2,croEbene.y2);
    cr->stroke();

    cr->set_source_rgb(1,1,1);
    cr->arc(croEbene.x1,croEbene.y1,5,0,2*M_PI);
    cr->fill();
    cr->arc(croEbene.x2,croEbene.y2,5,0,2*M_PI);
    cr->fill();

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

    cr->set_source_rgb(0,0,0);
    draw_text(cr, croGP.P123().x,croGP.P123().y+27, sId+"123");

    cr->set_source_rgb(1,1,1);
    cr->arc(G0.x,G0.y,25,0,2*M_PI); cr->fill();
    cr->set_source_rgb(0,0,0);
    cr->arc(G0.x,G0.y,25,0,2*M_PI); cr->stroke();

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

    for ( int i{0}; i<3; ++i)
	{
	SPoint const GPoint = croGP.GPoint(i);

    	cr->set_line_cap(Cairo::LINE_CAP_ROUND);
    	cr->set_line_width(4);
    	cr->set_source_rgba(.5,.5,.5,.5);
    	cr->move_to(G0.x,G0.y);
    	cr->line_to(GPoint.x,GPoint.y);

    	if (nId == 2)
    	    {
    	    cr->line_to(g_vGrundPunkte[0].GPoint(i).x,g_vGrundPunkte[0].GPoint(i).y);
    	    }
    	cr->stroke();

    	if (i==0)
    	    {
	    cr->set_source_rgb(0,0,0);
	    draw_text(cr, (G0.x+GPoint.x)/2,(G0.y+GPoint.y)/2,
          		   std::to_string((int) CalcDistance(G0, GPoint) ));

	    cr->set_source_rgba(.5,.5,.5,.5);
    	    }
	}
//--------------------------------
    for ( int i{0}; i<3; ++i)
	{
	SPoint const GPoint = croGP.GPoint(i);

	cr->set_source_rgb(.8,.8,1);
	cr->arc(GPoint.x,GPoint.y,10,0,2*M_PI);
	cr->fill();
	cr->set_source_rgb(0,0,0);
	cr->set_line_width(2);
	cr->arc(GPoint.x,GPoint.y,10,0,2*M_PI);
	cr->stroke();

    	cr->set_source_rgb(0,0,0);
    	draw_text(cr, GPoint.x,GPoint.y, sId+std::to_string(i+1));
	}

    cr->set_source_rgb(0,0,0);
    draw_text(cr, G0.x,G0.y-42, sId+"0");
    }

void draw_dreieck(Cairo::RefPtr<Cairo::Context> const & cr,
	          VPolDreieck const & croPD, int nId)
    {
    cr->set_line_cap(Cairo::LINE_CAP_ROUND);
    cr->set_line_width(13);
    cr->set_source_rgba(0,1,0,.5);
    cr->move_to(croPD[0].x,croPD[0].y);
    cr->line_to(croPD[1].x,croPD[1].y);
    cr->line_to(croPD[2].x,croPD[2].y);
    cr->fill();
    cr->set_source_rgb(0,0,0);
    cr->set_line_width(2);
    cr->move_to(croPD[0].x,croPD[0].y);
    cr->line_to(croPD[1].x,croPD[1].y);
    cr->line_to(croPD[2].x,croPD[2].y);
    cr->line_to(croPD[0].x,croPD[0].y);
    cr->stroke();

    SUmkreis tUmkreis = Umkreis( croPD[0], croPD[1], croPD[2] );
    cr->set_source_rgba(1,1,1,.25);
    cr->arc(tUmkreis.MidPnt.x,tUmkreis.MidPnt.y,tUmkreis.Radius,0,2*M_PI);
    cr->fill();
    cr->set_source_rgb(0,0,0);
    cr->arc(tUmkreis.MidPnt.x,tUmkreis.MidPnt.y,tUmkreis.Radius,0,2*M_PI);
    cr->stroke();

    cr->set_source_rgb(0,0,0);
    draw_text(cr,  croPD[0].x,croPD[0].y, "P12");
    draw_text(cr,  croPD[1].x,croPD[1].y, "P13");
    draw_text(cr,  croPD[2].x,croPD[2].y, "P23");
    draw_text(cr,  croPD[0].x,croPD[0].y, "");
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
	w0 = dTransX = width/2;
	h0 = dTransY = height/2;
	bTransInitialized = true;
	}
    if ( (w0!=width/2) || (h0!=height/2) )
	{
	dTransX -= w0 - width/2;  w0 = width/2;
	dTransY -= h0 - height/2; h0 = height/2;
	}

//a    cr->save();

    Cairo::Matrix matrix(1.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    matrix.scale(dScale,dScale);
    matrix.translate(dTransX/dScale, dTransY/dScale);
    cr->transform(matrix);
/*
    cr->set_source_rgb(.8, .8, 1);
    cr->rectangle(-150, -150, 300, 300);
    cr->fill();
*/

// static content
    if (g_vEbenenLagen.size() == 3)
	{
//	cr->set_source_rgba(.0, 1.0, .0, .25);
	g_vPolDreieck[0] = ( CalcPolpunkt(g_vEbenenLagen[1-1], g_vEbenenLagen[2-1]) );
	g_vPolDreieck[1] = ( CalcPolpunkt(g_vEbenenLagen[1-1], g_vEbenenLagen[3-1]) );
	g_vPolDreieck[2] = ( CalcPolpunkt(g_vEbenenLagen[2-1], g_vEbenenLagen[3-1]) );
	draw_dreieck(cr, g_vPolDreieck, 0);
	}

// interactive content
    if ( g_tCollision.eWhat == SCollision::EWhat::Ebene )
	{
	cr->set_source_rgb(.5,.5,0);
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
//a    cr->restore();

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

	auto const S = sin(t*2*M_PI);
	auto const C = cos(t*2*M_PI);
	auto const & A0{g_vGrundPunkte[0].G0()};
	auto const & A1 {g_vGrundPunkte[0].GPoint(0)};
	auto const & A2 {g_vGrundPunkte[0].GPoint(1)};
	auto const & A3 {g_vGrundPunkte[0].GPoint(2)};
	auto const & B1 {g_vGrundPunkte[1].GPoint(0)};
	auto const & B2 {g_vGrundPunkte[1].GPoint(1)};
	auto const & B3 {g_vGrundPunkte[1].GPoint(2)};
	auto const & B0{g_vGrundPunkte[1].G0()};

	auto const GL = CalcDistance( A0, B0 );
	auto const AL = CalcDistance( A0, A1 );
	auto const BL = CalcDistance( B0, B1 );
	auto const CL = CalcDistance( A1, B1 );
	auto const DL = CalcDistance( A1, B0 );

	auto const pa = SPoint{A0.x+AL*S, A0.y+AL*C};
	auto const  d = CalcDistance( pa , B0 );

	auto const epsi = M_PI-CalcAlpha((pa.x-B0.x), -(pa.y-B0.y), d);
	auto const beta = CalcAlpha(BL,CL,d);
	auto const gama = CalcAlpha(CL,BL,d);

	SPoint pb  = SPoint{B0.x+BL*sin(epsi+gama), B0.y+BL*cos(epsi+gama)};
	SPoint pab;
	if ( m_bDurchschlagen )
	    {
	    pb  = SPoint{B0.x+BL*sin(epsi-gama), B0.y+BL*cos(epsi-gama)};
	    pab = SPoint{pa.x+CL*sin(epsi+beta+M_PI), pa.y+CL*cos(epsi+beta+M_PI)};
	    }
	else
	    {
	    pb  = SPoint{B0.x+BL*sin(epsi+gama), B0.y+BL*cos(epsi+gama)};
	    pab = SPoint{pa.x+CL*sin(epsi-beta+M_PI), pa.y+CL*cos(epsi-beta+M_PI)};
	    }
//	std::cout << epsi/M_PI*180 << ", a="  << (B0.y-pa.y) << ", b="  << (B0.x-pa.x) << ", c="  << d << ", " << '\n';

	cr->set_line_width(lw+2);
	cr->set_source_rgb(1,1,1);
	if ( d > (BL+CL) )
	    cr->set_source_rgb(0,0,0);

	cr->move_to(A0.x, A0.y);
	cr->line_to(pa.x, pa.y);
	cr->stroke();
/*
	cr->set_source_rgb(.4, .4, 1);
	cr->move_to(B0.x, B0.y);
	cr->line_to(B0.x+d*sin(epsi), B0.y+d*cos(epsi));
	cr->stroke();
*/
	cr->move_to(pa.x, pa.y);
//	cr->line_to(pa.x+CL*sin(epsi-beta+M_PI), pa.y+CL*cos(epsi-beta+M_PI));
	cr->line_to(pab.x, pab.y);
	cr->stroke();

	cr->move_to(B0.x, B0.y);
	cr->line_to(pb.x, pb.y);
	cr->stroke();
/*
	Durchschlag
	cr->set_source_rgb(0.5,0.5,0.5);
	cr->move_to(B0.x, B0.y);
	cr->line_to(pb2.x, pb2.y);
	cr->stroke();

	cr->move_to(pa.x, pa.y);
	cr->line_to(pa.x+CL*sin(epsi+beta+M_PI), pa.y+CL*cos(epsi+beta+M_PI));
	cr->stroke();
*/

	if ( CalcDistance(A1, pa)<20 )
	    if ( CalcDistance(B1, pb)<30 )
		cr->set_source_rgb(.7, 1, .7);
	    else
		cr->set_source_rgb(1, .4, .7);
	else if ( CalcDistance(A2, pa)<20 )
	    if ( CalcDistance(B2, pb)<30 )
		cr->set_source_rgb(.7, 1, .7);
	    else
		cr->set_source_rgb(1, .4, .7);
	else if ( CalcDistance(A3, pa)<20 )
	    if ( CalcDistance(B3, pb)<30 )
		cr->set_source_rgb(.7, 1, .7);
	    else
		cr->set_source_rgb(1, .4, .7);
	else
	    cr->set_source_rgb(0,0,0);

	if ( d > (BL+CL) )
	    cr->set_source_rgb(1,0,0);

	cr->set_line_width(lw);

	cr->move_to(A0.x, A0.y);
	cr->line_to(pa.x, pa.y);
	cr->stroke();
/*
	cr->set_source_rgb(.4, .4, 1);
	cr->move_to(B0.x, B0.y);
	cr->line_to(B0.x+d*sin(epsi), B0.y+d*cos(epsi));
	cr->stroke();
*/
	cr->move_to(pa.x, pa.y);
//	cr->line_to(pa.x+CL*sin(epsi-beta+M_PI), pa.y+CL*cos(epsi-beta+M_PI));
	cr->line_to(pab.x, pab.y);
	cr->stroke();

	cr->move_to(B0.x, B0.y);
	cr->line_to(pb.x, pb.y);
	cr->stroke();

	cr->arc(pa.x,pa.y,3,0,2*M_PI);
	cr->stroke();

	cr->arc(pb.x,pb.y,3,0,2*M_PI);
	cr->stroke();

	cr->set_source_rgb(1,1,1);

	cr->arc(pa.x,pa.y,3,0,2*M_PI);
	cr->fill();

	cr->arc(pb.x,pb.y,3,0,2*M_PI);
	cr->fill();



	// Ebene **************************************
	SPoint const ep1{g_vEbenenLagen[0].x1,g_vEbenenLagen[0].y1};
	SPoint const ep2{g_vEbenenLagen[0].x2,g_vEbenenLagen[0].y2};
	auto       we  = CalcVectorSlope(ep1, ep2);
	auto       le  = CalcDistance(ep1, ep2);

/*
 	auto       wg  = CalcVectorSlope(pa, pb);
	auto       wb  = CalcVectorSlope(pb, ep1);
	auto       dw  = CalcVectorDiff(pa, pb, ep1);
*/
	auto       wg  = CalcVectorSlope(pa, pb);
	auto       wb  = CalcVectorSlope(B1, ep1);
	auto       dw  = CalcVectorDiff(B1, A1, ep1);

	auto      dwe  = CalcVectorDiff(ep1, B1, ep2);


	auto const lve = CalcDistance(B1, ep1);

//	std::cout << " dw: " << dw/M_PI*180 << ",  wg: " << wg/M_PI*180 << ",  wb: " << wb/M_PI*180 << ",  lve: " << lve << '\n';

	SPoint be{pb.x+lve*sin(wg+dw-M_PI/2), pb.y+lve*cos(wg+dw-M_PI/2)};
	cr->set_source_rgb(0,1,1);
	cr->move_to(pb.x, pb.y);
	cr->line_to(be.x, be.y);
	cr->stroke();

	cr->set_source_rgb(1,0,0);
	cr->set_line_width(13);
	cr->move_to(be.x, be.y);
	cr->line_to(be.x+le*sin(wg+dw+dwe+M_PI/2), be.y+le*cos(wg+dw+dwe+M_PI/2));
	cr->stroke();

/*
	cr->move_to(pa.x, pa.y);
	cr->line_to(ep1.x, ep1.y);
	cr->line_to(pb.x, pb.y);
	cr->line_to(pa.x, pa.y);
	cr->fill();
*/
	}





    //  GUI
    auto constexpr uiOffset{ 32.0};
    auto constexpr uiBaseWd{  8.0};
    auto constexpr uiBaseLn{128.0};
/*
    g_dAnimateMax;
    g_dAnimate;
    g_dAnimateMin;
*/
    auto const dAniGui = (uiBaseLn-2*uiBaseWd)*g_dAnimate/(g_dAnimateMax-g_dAnimateMin);

/*
    cr->set_source_rgb(.8, .8, .9);
    cr->rectangle((uiOffset/2-dTransX)/dScale, (uiOffset/2-dTransY)/dScale, (uiOffset/2*3)/dScale, (uiBaseLn+uiOffset)/dScale);
    cr->fill();
*/

    // Buttons
    int i{0};
    for ( auto const & a:m_voButtons )
	{
	++i;
	cr->set_source_rgb(.8, .8, .9);
	if ( a.Collision({gx,gy}) )
	    cr->set_source_rgb(0,1,0);
	cr->rectangle( (a.x-dTransX)/dScale, (a.y-dTransY)/dScale, (a.w)/dScale, (a.h)/dScale);
	cr->fill();

	cr->set_source_rgb(0,0,0);
	draw_text(cr,  (a.x+a.w/2-dTransX)/dScale, (a.y+a.h/2-dTransY)/dScale, a.t, dScale);
	}



    cr->set_line_cap(Cairo::LINE_CAP_ROUND);

    cr->set_source_rgb(0,0,0);

    cr->set_line_width(uiBaseWd/dScale);

    cr->move_to((uiOffset-dTransX)/dScale, (uiOffset-dTransY         )/dScale);
    cr->line_to((uiOffset-dTransX)/dScale, (uiOffset-dTransY+uiBaseLn)/dScale);
    cr->stroke();


    cr->set_line_width(3.0/dScale);

    cr->move_to((uiOffset-dTransX           )/dScale, (uiOffset-dTransY+uiBaseWd/2)/dScale);
    cr->line_to((uiOffset-dTransX+2*uiBaseWd)/dScale, (uiOffset-dTransY+uiBaseWd/2)/dScale);
    cr->stroke();

    cr->move_to((uiOffset-dTransX           )/dScale, (uiOffset-dTransY-uiBaseWd/2+uiBaseLn)/dScale);
    cr->line_to((uiOffset-dTransX+2*uiBaseWd)/dScale, (uiOffset-dTransY-uiBaseWd/2+uiBaseLn)/dScale);
    cr->stroke();

    cr->set_source_rgb(0,1,0);
    cr->set_line_width(1.5/dScale);

    cr->move_to((uiOffset-dTransX           )/dScale, (uiOffset-dTransY-uiBaseWd/2+uiBaseLn-dAniGui)/dScale);
    cr->line_to((uiOffset-dTransX+2*uiBaseWd)/dScale, (uiOffset-dTransY-uiBaseWd/2+uiBaseLn-dAniGui)/dScale);
    cr->stroke();

/*
    cr->set_source_rgb(1,1,0);
    cr->arc(gx,gy,3,0,2*M_PI);
    cr->fill();
*/
    return true;
    }
