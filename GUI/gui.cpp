//#include <gtkmm/button.h>
#include <gtkmm/hvscale.h>
#include <gtkmm/main.h>
#include <gtkmm/stock.h>
#include <gtkmm/window.h>
#include <sys/stat.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    Gtk::Main app(argc, argv);
    Gtk::Window window;

    window.set_title("Texture Synthesis");
	struct stat buffer;
	if(stat("ims/1.png", &buffer) == 0)
        window.set_icon_from_file("ims/1.png");
    window.maximize();

    // Gtk::Button but(Gtk::Stock::QUIT);
    // window.add(but);
    // but.show();
    // but.signal_clicked().connect([]() { Gtk::Main::quit(); });

    Gtk::HScale slider(0.0, 1.0, 0.01);
    window.add(slider);
    slider.show();
    slider.signal_value_changed().connect([&slider]() { cout << slider.get_value() << endl; });

    Gtk::Main::run(window);
    return 0;
}
