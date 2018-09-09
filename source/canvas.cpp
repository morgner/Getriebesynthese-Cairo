#include "canvas.h"
#include "grundpunkt.h"
#include <string>
#include <array>

VEbenenLagen    g_vEbenenLagen;// E1,  E2,  E3  = 3 homologe Lagen einer Ebene
VPolDreieck     g_vPolDreieck; // P12, P13, P23 = 3 Polpunkte zu den o.g.Ebenenlagen
SUmkreis       g_tUmkreis;
VGelenke        g_vGelenke;
SAllCollision   g_tCollision{SAllCollision::EWhat::none, 0, 0};
CGrundpunkt     g_oGrundpunkt({-1, -1});
VGrundpunkte    g_vGrundPunkte;

double dScale  {1.0}; // scale
bool   ti {false}; // t-initialized?
double dTransX {0.0}; // translation x
double dTransY {0.0}; // translation y
double dTransXStart{0.0}; // translation x
double dTransYStart{0.0}; // translation y
double cx {0.0}; // mouse x
double cy {0.0}; // mouse y

double dColliderMouseDistance{0.0}; // collider distance

double dcdx{0}; // distance to a specified point
double dcdy{0};

double dmosx{0}; // mouse offset scaled
double dmosy{0};

template<typename A, typename B>
    auto CalcDistance(A const & a, B const & b)
	{
	return sqrt( pow((a.x-b.x),2) + pow((a.y-b.y),2) );
	}

double MouseDistance( SPoint const & p, SPoint const & pMouse )
    {
    dcdx = (p.x - pMouse.x)*dScale;
    dcdy = (p.y - pMouse.y)*dScale;
    return CalcDistance(p, pMouse);
    }

SAllCollision MouseCollision(SPoint const & crtMousePoint)
    {
    SAllCollision ac;
    ac.eWhat  = SAllCollision::EWhat::none;
    int cnt{0};

    dColliderMouseDistance = 12.0;
    for (auto const & a:g_vEbenenLagen)
	{
	if ( MouseDistance({a.x1, a.y1}, crtMousePoint) < dColliderMouseDistance )
	    {
	    ac.eWhat  = SAllCollision::EWhat::Ebene;
	    ac.nIndex = cnt;
	    ac.nSubIx = 0;
	    dmosx = dcdx;
	    dmosy = dcdy;
	    break;
	    }
	if ( MouseDistance({a.x2, a.y2}, crtMousePoint) < dColliderMouseDistance )
	    {
	    ac.eWhat  = SAllCollision::EWhat::Ebene;
	    ac.nIndex = cnt;
	    ac.nSubIx = 1;
	    dmosx = dcdx;
	    dmosy = dcdy;
	    break;
	    }
	SPoint m = a.M();
	if ( MouseDistance({m.x, m.y}, crtMousePoint) < dColliderMouseDistance )
	    {
	    ac.eWhat  = SAllCollision::EWhat::Ebene;
	    ac.nIndex = cnt;
	    ac.nSubIx = 2;
	    dmosx = dcdx;
	    dmosy = dcdy;
	    break;
	    }
	++cnt;
	}
    cnt = 0;
    for (auto const & a:g_vGrundPunkte)
	{
	if ( MouseDistance( a.P123(), crtMousePoint) < dColliderMouseDistance )
	    {
	    ac.eWhat  = SAllCollision::EWhat::Grundpunkt;
	    ac.nIndex = cnt;
	    ac.nSubIx = 0;
	    dmosx = dcdx;
	    dmosy = dcdy;
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


CCanvas::CCanvas()
    : m_bFirstClick(false),
      x1(0),y1(0),
      x2(0),y2(0)
    {
    add_events(Gdk::BUTTON_PRESS_MASK | Gdk::SCROLL_MASK | Gdk::SMOOTH_SCROLL_MASK);
    add_events(Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK);
    add_events(Gdk::BUTTON1_MOTION_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::POINTER_MOTION_MASK);
    add_events(Gdk::KEY_PRESS_MASK | Gdk::KEY_RELEASE_MASK);
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
    auto const ex = event->x/dScale-dTransX/dScale;
               cx = event->x;
               dTransXStart = dTransX;
    auto const ey = event->y/dScale-dTransY/dScale;
               cy = event->y;
               dTransYStart = dTransY;

    if( (event->type == GDK_BUTTON_PRESS) && (event->button == 1) )
	{
        if ( !m_bFirstClick )
            {
            x1=x2=ex;
            y1=y2=ey;
            }
        if (g_tCollision.eWhat == SAllCollision::EWhat::none)
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

bool CCanvas::on_motion_notify_event(GdkEventMotion *event)
    {
    auto const ex=event->x/dScale-dTransX/dScale;
    auto const ey=event->y/dScale-dTransY/dScale;

    if ( event->type & GDK_MOTION_NOTIFY )
	if ( event->state & GDK_BUTTON1_MASK )
	    {
	    switch (g_tCollision.eWhat)
		case SAllCollision::EWhat::none:
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
			    dTransX = dTransXStart - (cx-event->x);
			    dTransY = dTransYStart - (cy-event->y);
			    break;
			}
		    break;

		case SAllCollision::EWhat::Ebene:
		    MoveEbenenPunkt(ex+dmosx/dScale, ey+dmosy/dScale, m_nLenEbene);
		    break;

		case SAllCollision::EWhat::Grundpunkt:
		    g_vGrundPunkte[g_tCollision.nIndex].Update({ex+dcdx/dScale, ey+dcdy/dScale});
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
    if (g_tCollision.eWhat == SAllCollision::EWhat::none)
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
    dScale *= (event->delta_y<0)?.9:1.1; if (dScale<.01) dScale=.01;
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
	       std::string const & crsText)
{
  cr->save();

  // http://developer.gnome.org/pangomm/unstable/classPango_1_1FontDescription.html
  Pango::FontDescription font;

  font.set_family("Monospace");
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
    cr->save();
    cr->set_line_cap(Cairo::LINE_CAP_ROUND);
    cr->set_line_width(13);
    cr->set_source_rgb(1,0,0);
    cr->move_to(croEbene.x1,croEbene.y1);
    cr->line_to(croEbene.x2,croEbene.y2);
    cr->stroke();
    cr->restore();

    cr->save();
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
    cr->restore();
    }

void draw_grundpunkt(Cairo::RefPtr<Cairo::Context> const & cr,
	             CGrundpunkt & croGP, int nId)
    {
    croGP.Show( cr, g_vPolDreieck);

    SPoint const G0 = croGP.G0();

    std::string sId = std::array<std::string,5>{"U","A","B","C","E"}[(nId>4)?4:nId];

    cr->save();
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
    cr->restore();

    for ( int i{0}; i<3; ++i)
	{
	SPoint const GPoint = croGP.GPoint(i);

    	cr->save();
    	cr->set_line_cap(Cairo::LINE_CAP_ROUND);
    	cr->set_line_width( (i==0)?16:4);
    	cr->set_source_rgba(.5,.5,.5,.5);
    	cr->move_to(G0.x,G0.y);
    	cr->line_to(GPoint.x,GPoint.y);

    	if (nId == 2)
    	    {
    	    cr->line_to(g_vGrundPunkte[0].GPoint(i).x,g_vGrundPunkte[0].GPoint(i).y);
    	    }
    	cr->stroke();
    	cr->restore();

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

	cr->save();
	cr->set_source_rgb(.8,.8,1);
	cr->arc(GPoint.x,GPoint.y,10,0,2*M_PI);
	cr->fill();
	cr->set_source_rgb(0,0,0);
	cr->set_line_width(2);
	cr->arc(GPoint.x,GPoint.y,10,0,2*M_PI);
	cr->stroke();
	cr->restore();

    	cr->save();
    	cr->set_source_rgb(0,0,0);
    	draw_text(cr, GPoint.x,GPoint.y, sId+std::to_string(i+1));
    	cr->restore();
	}
//--------------------------------

    cr->set_source_rgb(0,0,0);
    draw_text(cr, G0.x,G0.y-42, sId+"0");

    }

void draw_dreieck(Cairo::RefPtr<Cairo::Context> const & cr,
	          VPolDreieck const & croPD, int nId)
    {
    cr->save();
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
    cr->restore();

    SUmkreis tUmkreis = Umkreis( croPD[0], croPD[1], croPD[2] );
    cr->save();
    cr->set_source_rgba(1,1,1,.25);
    cr->arc(tUmkreis.MidPnt.x,tUmkreis.MidPnt.y,tUmkreis.Radius,0,2*M_PI);
    cr->fill();
    cr->set_source_rgb(0,0,0);
    cr->arc(tUmkreis.MidPnt.x,tUmkreis.MidPnt.y,tUmkreis.Radius,0,2*M_PI);
    cr->stroke();
    cr->restore();

    cr->save();
    cr->set_source_rgb(0,0,0);
    draw_text(cr,  croPD[0].x,croPD[0].y, "P12");
    draw_text(cr,  croPD[1].x,croPD[1].y, "P13");
    draw_text(cr,  croPD[2].x,croPD[2].y, "P23");
    draw_text(cr,  croPD[0].x,croPD[0].y, "");
    cr->restore();
    }

bool CCanvas::on_draw(Cairo::RefPtr<Cairo::Context> const & cr)
    {
    Gtk::Allocation allocation = get_allocation();
    const int width  = allocation.get_width();
    const int height = allocation.get_height();

    static auto w0 = width/2;
    static auto h0 = height/2;

	if ( false == ti )
	{
	w0 = dTransX = width/2;
	h0 = dTransY = height/2;
	ti = true;
	}
    if ( (w0!=width/2) || (h0!=height/2) )
	{
	dTransX -= w0 - width/2;  w0 = width/2;
	dTransY -= h0 - height/2; h0 = height/2;
	}

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
    if ( g_tCollision.eWhat == SAllCollision::EWhat::Ebene )
	{
	cr->save();
	cr->set_source_rgb(.5,.5,0);
	cr->set_line_width(11);
	switch ( g_tCollision.nSubIx )
	    {
	    case 0:
		cr->arc(g_vEbenenLagen[g_tCollision.nIndex].x1,
			g_vEbenenLagen[g_tCollision.nIndex].y1,dColliderMouseDistance,0,2*M_PI);
		break;
	    case 1:
		cr->arc(g_vEbenenLagen[g_tCollision.nIndex].x2,
			g_vEbenenLagen[g_tCollision.nIndex].y2,dColliderMouseDistance,0,2*M_PI);
		break;
	    case 2:
		cr->set_line_cap(Cairo::LINE_CAP_ROUND);
		cr->set_line_width(dColliderMouseDistance*2);
		cr->move_to(g_vEbenenLagen[g_tCollision.nIndex].x1,g_vEbenenLagen[g_tCollision.nIndex].y1);
		cr->line_to(g_vEbenenLagen[g_tCollision.nIndex].x2,g_vEbenenLagen[g_tCollision.nIndex].y2);
		cr->stroke();
		break;
	    }
	cr->fill();
	cr->restore();
	}

    if ( g_tCollision.eWhat == SAllCollision::EWhat::Grundpunkt )
	{
	cr->save();
	cr->set_source_rgb(.5,.5,0);
	cr->set_line_width(11);
	cr->arc(g_vGrundPunkte[g_tCollision.nIndex].P123().x,
		g_vGrundPunkte[g_tCollision.nIndex].P123().y, 21, 0, 2*M_PI);
	cr->fill();
	cr->restore();
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
    cr->save();
    cr->set_source_rgb(.8,.8,.0);
    cr->set_line_width(2);
    cr->arc(mx/s-tx/s, my/s-ty/s, 7/s,0,2*M_PI);
    cr->fill();
    cr->restore();

    cr->save();
    cr->set_source_rgb(.0,.8,.8);
    cr->set_line_width(2);
    cr->arc(-tx/s, -ty/s, 11/s,0,2*M_PI);
    cr->fill();
    cr->restore();

    cr->set_source_rgb(0,0,0);
    cr->set_font_size(24.0/s);
    draw_text(cr, 200/s-tx/s, 50/s-ty/s, std::to_string((int)tx)   + ", " + std::to_string((int)ty));
    draw_text(cr, 200/s-tx/s, 65/s-ty/s, std::to_string((int)(tx/s)) + ", " + std::to_string((int)(ty/s)));
    draw_text(cr, 200/s-tx/s, 80/s-ty/s, std::to_string((int)mx)   + ", " + std::to_string((int)my));
    draw_text(cr, 200/s-tx/s, 95/s-ty/s, std::to_string((int)(mx/s)) + ", " + std::to_string((int)(my/s)));
*/
    return true;
    }

