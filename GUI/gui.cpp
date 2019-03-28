//#include <gtkmm/button.h>
// #include <gtkmm/alignment.h>
// #include <gtkmm/hvscale.h>
// #include <gtkmm/hvbuttonbox.h>
// #include <gtkmm/main.h>
// #include <gtkmm/stock.h>
// #include <gtkmm/window.h>
// #include <gtkmm/action.h>
// #include <gtkmm/hvbox.h>
// #include <gtkmm/range.h>
#include <gtkmm.h>

#include "../synthetizer.hpp"
#include "../stb_image_write.h"

#include <sys/stat.h>
#include <iostream>

using namespace std;

class ScalingImage : public Gtk::Image {
public:
	explicit ScalingImage(Glib::RefPtr<Gdk::Pixbuf> pixbuf, Gdk::InterpType interp = Gdk::INTERP_BILINEAR):
		Gtk::Image(pixbuf),
		m_original(pixbuf),
		m_interp(interp) {}

	void setPixbuf(Glib::RefPtr<Gdk::Pixbuf> pixbuf) {
		m_original = pixbuf;
		double scale =  min(get_pixbuf()->get_width() / double(m_original->get_width()),
							get_pixbuf()->get_height() / double(m_original->get_height()));
		Glib::RefPtr<Gdk::Pixbuf> scaled =
			m_original->scale_simple(max(1.0, scale*m_original->get_width()), max(1.0, scale*m_original->get_height()), m_interp);
		Gtk::Image::set(scaled);
	}

protected:
	virtual void on_size_allocate(Gtk::Allocation & r) {
		double scale =  min((r.get_width()-24) / double(m_original->get_width()),
							(r.get_height()-24) / double(m_original->get_height()));
		Glib::RefPtr<Gdk::Pixbuf> scaled =
			m_original->scale_simple(max(1.0, scale*m_original->get_width()), max(1.0, scale*m_original->get_height()), m_interp);
		Gtk::Image::set(scaled);
		Gtk::Image::on_size_allocate(r);
	}
  
private:
	Glib::RefPtr<Gdk::Pixbuf> m_original;
	Gdk::InterpType m_interp;

};

static ScalingImage *image = NULL;
static Gtk::HBox *hb;

// Parameters
static int W = 5, H = 3;
VD r = VD(9, 0);
static int c = 3;
static double kappa = 0.2;
static bool compute_co = false;
static bool to_tor = false;
static char const *filename = "ims/2.png";
static uchar *E, *Ed, *E2;
static int m, md, m2;
static double is_tore, new_E;
static bool have_folder;
static char folder[100];
static Pix *Ss[9];
static int Ws[9], Hs[9];
static uchar *El[9];
static int l = 0, L;

void clean() {
	if(new_E) delete[] E;
	else free(E);
	if(md < m) delete[] Ed;
	if(!is_tore) delete[] E2;
	for(int i = 0; i < 9; i++) {
		if(!Ss[i]) delete[] Ss[i];
		if(!El[i]) delete[] El[i];
	}
}

void update() {
	if(l > L) {
		if(md < m) {
			int Wh, Hh;
			uchar* Sh = magnify(md, E, m, Ss[L], W, H, Wh, Hh);
			stbi_write_png("out.png", Wh, Hh, 3, Sh, 0);
			delete[] Sh;
		} else
			save_smooth(Ss[L], Ws[L], Hs[L], E, m, "out.png");
		image->setPixbuf(Gdk::Pixbuf::create_from_file("out.png"));
		l = -1;
		return;
	}
	l ++;
	synthesize_step(l-1, Ss, Ws, Hs, E2, El, m, m2,
					r, L, have_folder, folder, compute_co, c, kappa);
	if(image != NULL)
		image->setPixbuf(Gdk::Pixbuf::create_from_file("out.png"));
	else {
		image = new ScalingImage(Gdk::Pixbuf::create_from_file("out.png"));
		hb->pack_start(*image);
		image->show();
	}
	update();
}

void slider_fun(Gtk::HScale *slider, int i) {
	r[i] = slider->get_value();
	if(l == -1) {
		l = i;
		update();
	} else if(i < l) l = i;
}

int main(int argc, char* argv[]) {
	Gtk::Main app(argc, argv);
	Gtk::Window window;

	window.set_title("Texture Synthesis");
	struct stat buffer;
	if(stat("ims/1.png", &buffer) == 0)
		window.set_icon_from_file("ims/1.png");
	window.set_default_size(600, 400);

	hb = new Gtk::HBox(false, 0);
	window.add(*hb);
	Gtk::VButtonBox box(Gtk::BUTTONBOX_SPREAD);
	hb->pack_end(box);

	Gtk::HScale sliders[9];
	for(int i = 0; i < 9; i++) {
		sliders[i] = Gtk::HScale(0.0, 1.01, 0.01);
		box.pack_start(sliders[i]);
		Gtk::HScale *slider = sliders+i;
		sliders[i].signal_button_release_event().connect([slider, i](GdkEventButton *e) {
			Glib::Thread::create([slider, i]() {
				slider_fun(slider, i);
			}, true);
			return false;
		});
	}

	window.show_all();

	// Initialisation
	if(load_image(filename, to_tor, E, m, is_tore, new_E, Ed, md))
		return 1;
	init_variables(Ed, md, W, H, !is_tore, filename,
					m2, E2, have_folder, folder, L);
	init_live(W, H, E2, md, m2, L,
				Ss, Ws, Hs, El);
	Glib::Thread::create([]() {	update(); }, true);
	
	Gtk::Main::run(window);
	clean();
	return 0;
}



	// Gtk::Button but(Gtk::Stock::QUIT);
	// window.add(but);
	// but.show();
	// but.signal_clicked().connect([]() { Gtk::Main::quit(); });

	// Gtk::Alignment rt_al(Gtk::ALIGN_END, Gtk::ALIGN_START, 0, 0);
	// window.add(rt_al);