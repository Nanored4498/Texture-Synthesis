#include <gtkmm/hvscale.h>
#include <gtkmm/hvbuttonbox.h>
#include <gtkmm/label.h>
#include <gtkmm/button.h>
#include <gtkmm/main.h>
#include <gtkmm/stock.h>
#include <gtkmm/menubar.h>
#include <gtkmm/menu.h>
#include <gtkmm/imagemenuitem.h>
#include <gtkmm/filechooserdialog.h>
#include <gtkmm/window.h>

#include <glibmm/thread.h>
#include <glibmm/fileutils.h>

#include "../synthetizer.hpp"
#include "../stb_image_write.h"

#include <iostream>
#include <sys/stat.h>

using namespace std;

static int screen_width, screen_height;

class ScalingImage : public Gtk::Image {
public:
	explicit ScalingImage(Glib::RefPtr<Gdk::Pixbuf> pixbuf, Gdk::InterpType interp = Gdk::INTERP_BILINEAR):
		Gtk::Image(pixbuf),
		m_original(pixbuf),
		m_interp(interp) {}

	void setPixbuf(Glib::RefPtr<Gdk::Pixbuf> pixbuf) {
		m_original = pixbuf;
		double scale =  min(min((double) get_pixbuf()->get_width(), 0.8*screen_width) / double(m_original->get_width()),
							min((double) get_pixbuf()->get_height(), 0.9*screen_height) / double(m_original->get_height()));
		Glib::RefPtr<Gdk::Pixbuf> scaled =
			m_original->scale_simple(max(1.0, scale*m_original->get_width()), max(1.0, scale*m_original->get_height()), m_interp);
		Gtk::Image::set(scaled);
	}

	void setPixbuf(const std::string &filename) {
		try {
			Glib::RefPtr<Gdk::Pixbuf> pixb = Gdk::Pixbuf::create_from_file(filename);
			this->setPixbuf(pixb);
		} catch(const Glib::FileError &ex) {
			std::cerr << "FileError: " << ex.what() << std::endl;
		} catch(const Gdk::PixbufError &ex) {
			std::cerr << "PixbufError: " << ex.what() << std::endl;
		}
	}

protected:
	virtual void on_size_allocate(Gtk::Allocation & r) {
		double scale =  min(min((double) r.get_width()-24, 0.8*screen_width) / double(m_original->get_width()),
							min((double) r.get_height()-24, 0.9*screen_width) / double(m_original->get_height()));
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
static int W = 2, H = 2;
// VD r = {0.2, 0.3, 0.35, 0.3, 0.2, 0.1, 0.1, 0.1, 0.0};
VD r = VD(9, 0);
static int c = 2;
static double kappa = 0.2;
static bool to_tor = false;
static uchar *E, *Ed, *E2;
static int m, md, m2;
static double is_tore, new_E;
static bool have_folder;
static char folder[100];
static Pix *Ss[9];
static int Ws[9], Hs[9];
static uchar *El[9];
static int l = -1, L;
static bool saveE = true;

void clean() {
	if(!E) {
		if(new_E) delete[] E;
		else free(E);
	}
	if(md < m && !Ed) delete[] Ed;
	if(!is_tore && !E2) delete[] E2;
	for(int i = 0; i < 9; i++) {
		if(!Ss[i]) delete[] Ss[i];
		if(!El[i]) delete[] El[i];
	}
}

void update() {
	if(l > L) {
		if(!saveE) {
			l = -1;
			return;
		}
		if(md < m) {
			int Wh, Hh;
			uchar* Sh = magnify(md, E, m, Ss[L], W, H, Wh, Hh);
			stbi_write_png("out.png", Wh, Hh, 3, Sh, 0);
			delete[] Sh;
		} else save_smooth(Ss[L], Ws[L], Hs[L], E, m, "out.png");
		image->setPixbuf("out.png");
		if(l > L) {
			l = -1;
			return;
		} else update();
	}
	l ++;
	synthesize_step(l-1, Ss, Ws, Hs, E2, El, md, m2,
					r, L, have_folder, folder, false, c, kappa, saveE);
	if(image != NULL) {
		image->setPixbuf("out.png");
	} else {
		image = new ScalingImage(Gdk::Pixbuf::create_from_file("out.png"));
		hb->pack_start(*image);
		image->show();
	}
	update();
}

void jitter_fun(Gtk::HScale *slider, int i) {
	r[i] = slider->get_value();
	if(l == -1) {
		l = i;
		update();
	} else if(i < l) l = i;
}

void dim_fun(int W0, int H0) {
	if(W0 > 0) W = W0;
	if(H0 > 0) H = H0;
	init_live_WH(L, W, H, Ss, Ws, Hs);
	if(l == -1) {
		l = 0;
		update();
	} else l = 0;
}

void c_fun(int c0) {
	c = c0;
	if(l == -1) {
		l = 0;
		update();
	} else l = 0;
}

void saveE_fun() {
	saveE = !saveE;
	if(l == -1) {
		l = L;
		update();
	} else l = std::min(l, L);
}

void open_image(const char *filename) {
	bool pred_saveE = saveE;
	while(l != -1) {
		l = 42;
		usleep(100000);
	}
	saveE = pred_saveE;
	clean();
	if(load_image(filename, to_tor, E, m, is_tore, new_E, Ed, md)) {
		std::cerr << "Impossible to load the chosen image !" << std::endl;
		return;
	}
	init_variables(Ed, md, !is_tore, filename,
					m2, E2, have_folder, folder, L);
	init_live(W, H, E2, m2, L,
				Ss, Ws, Hs, El);
	l = 0;
	Glib::Thread::create([]() {	update(); }, false);
}

int main(int argc, char* argv[]) {
	Gtk::Main app(argc, argv);
	Gtk::Window window;

	window.set_title("Texture Synthesis");
	struct stat buffer;
	if(stat("ims/1.png", &buffer) == 0)
		window.set_icon_from_file("ims/1.png");
	window.set_default_size(600, 400);
	Glib::RefPtr<Gdk::Screen> screen = Gdk::Screen::get_default();
	screen_width = screen->get_width();
	screen_height = screen->get_height();

	Gtk::VBox vb;
	window.add(vb);

	// Menu Bar
	Gtk::MenuBar bar;
	Gtk::MenuItem fileItem("_File", true);
	bar.append(fileItem);
	Gtk::Menu fileMenu;
	fileItem.set_submenu(fileMenu);	
	Gtk::ImageMenuItem openItem(Gtk::Stock::OPEN);
	openItem.signal_activate().connect([&window]() {
		Gtk::FileChooserDialog dialog(window, "Open a sample image");
		dialog.set_current_folder("./ims/");
		
		Glib::RefPtr<Gtk::FileFilter> imageFilter = Gtk::FileFilter::create();
		imageFilter->set_name("Image File");
		imageFilter->add_mime_type("image/*");
		dialog.add_filter(imageFilter);
		
		dialog.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
		dialog.add_button(Gtk::Stock::OPEN, Gtk::RESPONSE_OK);
		
		if(dialog.run() == Gtk::RESPONSE_OK)
			open_image(dialog.get_filename().c_str());
	});
	fileMenu.append(openItem);
	vb.pack_start(bar, Gtk::PACK_SHRINK);

	hb = new Gtk::HBox(false, 0);
	vb.add(*hb);

	Gtk::VButtonBox box(Gtk::BUTTONBOX_SPREAD);
	hb->pack_end(box);

	Gtk::Button typeBut("E/S");
	typeBut.signal_clicked().connect([]() { saveE_fun(); });
	box.pack_start(typeBut);

	Gtk::VBox Wbox(false, 0);
	Gtk::Label Wlab("W");
	Wbox.pack_start(Wlab);
	Gtk::HScale Wslider(1, 6, 1);
	Wslider.set_value(W);
	Wslider.signal_button_release_event().connect([&Wslider](GdkEventButton *e) {
		int W0 = Wslider.get_value();
		Glib::Thread::create([W0]() { dim_fun(W0, 0); }, false);
		return false;
	});
	Wbox.pack_start(Wslider);
	box.pack_start(Wbox);

	Gtk::VBox Hbox(false, 0);
	Gtk::Label Hlab("H");
	Hbox.pack_start(Hlab);
	Gtk::HScale Hslider(1, 6, 1);
	Hslider.set_value(H);
	Hslider.signal_button_release_event().connect([&Hslider](GdkEventButton *e) {
		int H0 = Hslider.get_value();
		Glib::Thread::create([H0]() { dim_fun(0, H0); }, false);
		return false;
	});
	Hbox.pack_start(Hslider);
	box.pack_start(Hbox);

	Gtk::VBox Cbox(false, 0);
	Gtk::Label Clab("Number of corrections");
	Cbox.pack_start(Clab);
	Gtk::HScale Cslider(0, 5, 1);
	Cslider.set_value(c);
	Cslider.signal_button_release_event().connect([&Cslider](GdkEventButton *e) {
		c_fun(Cslider.get_value());
		return false;
	});
	Cbox.pack_start(Cslider);
	box.pack_start(Cbox);

	Gtk::Label jitter_lab("Jitter");
	box.pack_start(jitter_lab);
	Gtk::HScale sliders[9];
	for(int i = 0; i < 9; i++) {
		sliders[i] = Gtk::HScale(0.0, 1.01, 0.01);
		sliders[i].set_value_pos(Gtk::POS_RIGHT);
		box.pack_start(sliders[i]);
		Gtk::HScale *slider = sliders+i;
		sliders[i].signal_button_release_event().connect([slider, i](GdkEventButton *e) {
			Glib::Thread::create([slider, i]() {
				jitter_fun(slider, i);
			}, false);
			return false;
		});
	}

	if(argc > 1)
		open_image(argv[1]);

	window.show_all();
	
	Gtk::Main::run(window);
	clean();
	return 0;
}