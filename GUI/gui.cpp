//#include <gtkmm/button.h>
#include <gtkmm/hvscale.h>
#include <gtkmm/hvbuttonbox.h>
#include <gtkmm/main.h>
#include <gtkmm/stock.h>
#include <gtkmm/window.h>
#include <gtkmm/action.h>
#include <gtkmm/hvbox.h>
// #include <gtkmm/alignment.h>

#include <sys/stat.h>
#include <iostream>

using namespace std;

struct ScalingImage : public Gtk::Image {

	explicit ScalingImage(Glib::RefPtr<Gdk::Pixbuf> pixbuf, Gdk::InterpType interp = Gdk::INTERP_BILINEAR):
		Gtk::Image(pixbuf),
		m_original(pixbuf),
		m_interp(interp) {}

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

int main(int argc, char* argv[]) {
	Gtk::Main app(argc, argv);
	Gtk::Window window;

	window.set_title("Texture Synthesis");
	struct stat buffer;
	if(stat("ims/1.png", &buffer) == 0)
		window.set_icon_from_file("ims/1.png");
	window.set_default_size(600, 400);

	// Gtk::Button but(Gtk::Stock::QUIT);
	// window.add(but);
	// but.show();
	// but.signal_clicked().connect([]() { Gtk::Main::quit(); });

	// Gtk::Alignment rt_al(Gtk::ALIGN_END, Gtk::ALIGN_START, 0, 0);
	// window.add(rt_al);

	Gtk::HBox hb(false, 0);
	window.add(hb);

	Gtk::VButtonBox box(Gtk::BUTTONBOX_SPREAD);
	hb.pack_end(box);

	ScalingImage apsi(Gdk::Pixbuf::create_from_file("ims/1/res.png"));
	hb.pack_start(apsi);

	Gtk::HScale sliders[8];
	for(int i = 0; i < 8; i++) {
		sliders[i] = Gtk::HScale(0.0, 1.0, 0.01);
		box.pack_start(sliders[i]);
		Gtk::HScale *slider = sliders+i;
		sliders[i].signal_value_changed().connect([slider, i]() {
			cout << i+1 << " " << slider->get_value() << endl;
			//im.set("ims/1.png");
		} );
	}

	window.show_all();
	Gtk::Main::run(window);
	return 0;
}
