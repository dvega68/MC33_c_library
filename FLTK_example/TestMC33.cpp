/*
	File: TestMC33.cpp
	Programmed by: David Vega - dvega@uc.edu.ve
	July 2019
	February 2020
	July 2020
	July 2021
	September 2023
	This is an open source code. The distribution and use rights are under the terms of the MIT license (https://opensource.org/licenses/MIT)
*/

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstring>

#if defined(_WIN32) && !defined(__CYGWIN__)
#include <io.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Radio_Round_Button.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Color_Chooser.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>
#include <FL/fl_ask.H>

/*
#include <marching_cubes_33.h> // put -lMC33 in the link libraries
/*/
#include "../source/marching_cubes_33.c"
//*/
#include "../source/MC33_util_grd.c"

void drawsurface(surface *S) {
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glVertexPointer(3, GL_FLOAT, 0, S->V);
	glNormalPointer(GL_FLOAT, 0, S->N);
	glColorPointer(3, GL_UNSIGNED_BYTE, 4, S->color);

	glDrawElements(GL_TRIANGLES, 3*S->nT, GL_UNSIGNED_INT, S->T);
	
	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
}

void drawdraftsurface(surface *S) {
	glDisable(GL_LIGHTING);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glVertexPointer(3, GL_FLOAT, 12*sizeof(float), S->V);
	glColorPointer(3, GL_UNSIGNED_BYTE, 16, S->color);

	glDrawArrays(GL_POINTS, 0, S->nV>>2);

	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
	glEnable(GL_LIGHTING);
}


using namespace std;

class gl_window : public Fl_Gl_Window {
	float bgr, bgg, bgb;
	float light_position[2][4];
	float MGL[16];
	float oldM[9], iniV[3];
	float Tx, Ty, scale_factor;
	int oldx, oldy;
	float scale_size;
	float x_c, y_c, z_c;
	GLenum face;
	surface *surf;
	int drawflag;

	void FixViewport(int w, int h) {
		glViewport(0, 0, w, h);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		if (w < h)
			glOrtho(-1.0f, 1.0f, -1.0f*h/w, 1.0f*h/w, -10.0f, 10.0f);
		else
			glOrtho(-1.0f*w/h, 1.0f*w/h, -1.0f, 1.0f, -10.0f, 10.0f);
		glMatrixMode(GL_MODELVIEW);
	}

	void set_glMaterial() {
		float mat_specular[4] = { 1.0f, 0.7f, 0.2f, 1.0f};
		float shininess = 50.0f;

		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, face == GL_FRONT_AND_BACK);
		glColorMaterial(face, GL_AMBIENT_AND_DIFFUSE);
		glMaterialfv(face, GL_SPECULAR, mat_specular);
		glMaterialf(face, GL_SHININESS, shininess);
	}

	void draw() {
		if (!valid()) {
			valid(1);
			FixViewport(w(), h());
			float light0_ambient[4] = {0.2f, 0.2f, 0.2f, 1.0f};
			float light0_diffuse[4] = {0.8f, 0.8f, 0.8f, 1.0f};

			float light1_ambient[4] = {0.1f, 0.1f, 0.1f, 1.0f};
			float light1_diffuse[4] = {0.6f, 0.6f, 0.6f, 1.0f};

		// light 0 settings
			glEnable(GL_LIGHT0);
			glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
			glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
			glLightfv(GL_LIGHT0, GL_SPECULAR, light0_diffuse);

		// light 1 settings
			glEnable(GL_LIGHT1);
			glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
			glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
			glLightfv(GL_LIGHT1, GL_SPECULAR, light1_diffuse);

			glEnable(GL_DEPTH_TEST);
			glEnable(GL_LIGHTING);
			glEnable(GL_COLOR_MATERIAL);
			set_glMaterial();
		}
		glLoadIdentity();
		glClearColor(bgr, bgg, bgb, 0.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glLightfv(GL_LIGHT0, GL_POSITION, light_position[0]);
		glLightfv(GL_LIGHT1, GL_POSITION, light_position[1]);
		MGL[12] = Tx*scale_size - x_c*MGL[0] - y_c*MGL[4] - z_c*MGL[8];
		MGL[13] = Ty*scale_size - x_c*MGL[1] - y_c*MGL[5] - z_c*MGL[9];
		MGL[14] = -x_c*MGL[2] - y_c*MGL[6] - z_c*MGL[10];
		glLoadMatrixf(MGL);

		if (surf) {
			if (drawflag)
				drawdraftsurface(surf);
			else
				drawsurface(surf);
		}
	}

	void mapToSphere(float x, float y, float *NewV) {
		float length = x*x + y*y;
		if (length > 1.0f) {
			float norm = 1.0f/sqrtf(length);
			NewV[0] = x * norm;
			NewV[1] = y * norm;
			NewV[2] = 0.0f;
		} else {
			NewV[0] = x;
			NewV[1] = y;
			NewV[2] = sqrtf(1.0f - length);
		}
	}

public:
	gl_window(int X,int Y,int W,int H, const char*L = 0) : Fl_Gl_Window(X, Y, W, H, L) {
		end();
		Tx = Ty = bgr = bgg = bgb = x_c = y_c = z_c = 0.0f;
		scale_factor = scale_size = 1.0f;
		set_light_pos(0, 1.5f, -0.5f, 2.0f);
		set_light_pos(1, -1.5f, 0.5f, 2.0f);
		reset_view();
	}

	void set_scale(float scale) {
		scale_size = scale;
	}

	void set_center(float x, float y, float z) {
		x_c = x; y_c = y; z_c = z;
	}

	void set_light_pos(int i, float x, float y, float z) {
		light_position[i][0] = x;
		light_position[i][1] = y;
		light_position[i][2] = z;
		light_position[i][3] = 0.0f;
	}

	void SetBackgroundColor(float r, float g, float b) {bgr = r; bgg = g; bgb = b; redraw();}
	void SetSurface(surface *s) {surf = s; redraw();}

	void reset_view() {
		float *f;
		for (f = MGL + 12; --f != MGL;) // MGL[0] is assigned later
			*f = 0.0f;
		for (f = oldM + 8; --f != oldM;) // oldM[0] and oldM[8] are assigned later
			*f = 0.0f;
		oldM[0] = oldM[4] = oldM[8] = MGL[0] = MGL[5] = MGL[10] = scale_factor = 1.0f;
		Tx = Ty = 0.0f;
		MGL[15] = scale_size;
		drawflag = 0;
		redraw();
	}

	void face_F_FB() {
		face = face == GL_FRONT? GL_FRONT_AND_BACK: GL_FRONT;
		set_glMaterial();
		redraw();
	}

	int handle(int event) {
		switch(event) {
			case FL_PUSH:
				oldx = Fl::event_x();
				oldy = Fl::event_y();
				if (Fl::event_button1()) {
						if (Fl::event_clicks()) {
							drawflag = 0;
							reset_view();
						}
						mapToSphere(2.0f*oldx/(w() - 1) - 1.0f,1.0f - 2.0f*oldy/(h() - 1), iniV);
				}
				return 1;
			case FL_DRAG:
			{
				int x = Fl::event_x();
				int y = Fl::event_y();
				switch(Fl::event_button()) {
					case FL_LEFT_MOUSE:
					{
						float qm[9]; // the 3 first elements of qm are the end vector to map to sphere
						// the elements 3 to 6 of qm are the rotation quaternion.
						// To simplify the calculation, the angle of rotation is multiplied by 2.
						mapToSphere(2.0f*x/(w() - 1) - 1.0f,1.0f - 2.0f*y/(h() - 1), qm);
						qm[4] = iniV[2]*qm[1] - iniV[1]*qm[2]; // qi
						qm[5] = iniV[0]*qm[2] - iniV[2]*qm[0]; // qj
						qm[6] = iniV[1]*qm[0] - iniV[0]*qm[1]; // qk
						float s = qm[4]*qm[4] + qm[5]*qm[5] + qm[6]*qm[6]; // sin^2(angle)
						if (s < 0.00001f*0.00001f) return 1;
						qm[3] = iniV[0]*qm[0] + iniV[1]*qm[1] + iniV[2]*qm[2]; // cos(angle) instead cos(angle/2)
						s += qm[3]*qm[3];
						s = (s > 0.0f ? 2.0f/s : 0.0f);
						// for notation see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
						float qi2 = qm[4]*s, qj2 = qm[5]*s, qk2 = qm[6]*s;
						float qiqr = qi2*qm[3], qjqr = qj2*qm[3], qkqr = qk2*qm[3];
						float qiqj = qm[4]*qj2, qiqk = qm[4]*qk2, qjqk = qm[5]*qk2;
						qi2 *= qm[4];
						qj2 *= qm[5];
						qk2 *= qm[6];
						// now qm is the rotation matrix:
						qm[0] = 1.0f - (qj2 + qk2); qm[1] = qiqj - qkqr; qm[2] = qiqk + qjqr;
						qm[3] = qiqj + qkqr; qm[4] = 1.0f - (qi2 + qk2); qm[5] = qjqk - qiqr;
						qm[6] = qiqk - qjqr; qm[7] = qjqk + qiqr; qm[8] = 1.0f - (qi2 + qj2);
						//Update the OpenGL transformation matrix:
						MGL[0] = qm[0]*oldM[0] + qm[3]*oldM[1] + qm[6]*oldM[2];
						MGL[4] = qm[0]*oldM[3] + qm[3]*oldM[4] + qm[6]*oldM[5];
						MGL[8] = qm[0]*oldM[6] + qm[3]*oldM[7] + qm[6]*oldM[8];
						MGL[1] = qm[1]*oldM[0] + qm[4]*oldM[1] + qm[7]*oldM[2];
						MGL[5] = qm[1]*oldM[3] + qm[4]*oldM[4] + qm[7]*oldM[5];
						MGL[9] = qm[1]*oldM[6] + qm[4]*oldM[7] + qm[7]*oldM[8];
						MGL[2] = qm[2]*oldM[0] + qm[5]*oldM[1] + qm[8]*oldM[2];
						MGL[6] = qm[2]*oldM[3] + qm[5]*oldM[4] + qm[8]*oldM[5];
						MGL[10] = qm[2]*oldM[6] + qm[5]*oldM[7] + qm[8]*oldM[8];
					}
					break;
					case FL_MIDDLE_MOUSE:
					{
						scale_factor += (y - oldy)*scale_factor*0.001f;
						MGL[15] = scale_factor*scale_size;
						oldy = y;
					}
					break;
					case FL_RIGHT_MOUSE:
					{
						float t = 2*(w() > h() ? scale_factor/h() : scale_factor/w());
						Tx += (x - oldx)*t;
						Ty -= (y - oldy)*t;
						oldx = x;
						oldy = y;
					}
					break;
				}
				drawflag = 1;
				redraw();
				return 1;
			}
			case FL_RELEASE:
				memcpy(oldM, MGL, 3*sizeof(float));
				memcpy(oldM + 3, MGL + 4, 3*sizeof(float));
				memcpy(oldM + 6, MGL + 8, 3*sizeof(float));
				if (drawflag) {
					if (--drawflag == 0)
						redraw();
				}
				return 1;
			case FL_MOUSEWHEEL:
				drawflag = 4;
				scale_factor += Fl::event_dy()*scale_factor*0.03125f;
				MGL[15] = scale_factor*scale_size;
				redraw();
				return 1;
			default:
				if (drawflag) {
					if (--drawflag == 0)
						redraw();
				}
				return(Fl_Gl_Window::handle(event));
		}
	}
};

int ask_replace_file(char *filename) {
	if (access(filename, F_OK))
		return 1;
	else
		return fl_choice("File already exists, overwrite?", "No", "Yes", 0);
}

long long filesize(const char *filename) {
	ifstream is (filename, ifstream::binary);
	is.seekg (0, is.end);
	return is.tellg();
}

class MyWindowClass : public Fl_Double_Window {
	gl_window *wgl;					// opengl window
	Fl_Box *numvert_box, *time_box, *info_box;
	Fl_Value_Input *isoval_input, *res_w[3], *ratiox, *ratioy, *ratioz;
	Fl_Double_Window *reswindow, *ratiowindow;
	int accept, defaultcolor;
	Fl_Check_Button *bigendian;
	vector <surface*> vs;
	size_t vs_i;
	_GRD *G;
	MC33 *M;

private:
	void fill_surface_info() {
		if (vs.empty()) {
			numvert_box->copy_label("");
			time_box->copy_label("");
			wgl->SetSurface(0);
		} else {
			ostringstream info;
			info << "Surface: " << vs_i + 1 << " | " << vs.size() << "   Vertices: ";
			info << vs[vs_i]->nV << "   Triangles: " << vs[vs_i]->nT;
			numvert_box->copy_label(info.str().c_str());
			info.str("");
			info << "Time: " << vs[vs_i]->user.i[0] << " ms";
			time_box->copy_label(info.str().c_str());
			isoval_input->value(vs[vs_i]->user.f[3]);
			wgl->SetSurface(vs[vs_i]);
		}
		numvert_box->redraw();
		time_box->redraw();
	}

	void CalcIso_cb2() {
		double iso = isoval_input->value();
		for (size_t i = 0; i != vs.size(); i++)
			if (vs[i]->user.f[3] == iso) {
				vs_i = i;
				fill_surface_info();
				return;
			}
		if (M) {
			int t = clock();
			surface *S = calculate_isosurface(M, iso);
			t = clock() - t;
			if (S) {
				if (S->nV > 0) {
					vs_i = vs.size();
					S->user.i[0] = t*1000/CLOCKS_PER_SEC;
					vs.push_back(S);
					fill_surface_info();
				} else {
					fl_message("Empty surface!");
					free_surface_memory(S);
				}
			} else
				fl_message("Insufficient memory!");
		} else
			fl_message("No grid data!");
	}

	static void CalcIso_cb(Fl_Widget *b, void *userdata) {
		MyWindowClass *w = reinterpret_cast<MyWindowClass*>(userdata);
		if (b == w->isoval_input && Fl::event() != FL_KEYBOARD) 
			return; // only enter key do callback in isoval_input
		w->CalcIso_cb2();
	}

	void Open_cb2() {
		int gridtype = 0;
		wgl->SetSurface(0);
		char *s = fl_file_chooser("Open grid file", "All files(*)", "");
		if (!s) {
			if (vs.size())
				wgl->SetSurface(vs[vs_i]);
			return;
		}
		free_memory_grd(G);
		G = 0;
		free_MC33(M);
		M = 0;
		char *c = strrchr(s, '.');
		if (c) {
			if (!strcmp(c, ".grb"))
				G = read_grd_binary(s);
			else if (!strcmp(c, ".grd"))
				G = read_grd(s);
			else if (!strcmp(c, ".dat")) {
				G = read_dat_file(s);
				gridtype = G != 0;
			}
		}
		if (!G) {
			long long fsize = filesize(s);
			gridtype = 1;
			switch (fsize) {
				case 8192:
					G = read_scanfiles(s, 64, 0);
					break;
				case 131072:
					G = read_scanfiles(s, 256, 1);
					break;
				case 524288:
					G = read_scanfiles(s, 512, 1);
					break;
				default:
				{
					unsigned int n[3];
					accept = 0;
					reswindow->show();
					while(reswindow->shown())
						Fl::wait();
					if (accept) {
						for(int i = 0; i != 3; ++i)
							n[i] = (int)res_w[i]->value();
						long long div = n[0]*n[1]*n[2];
						unsigned int nb = (div? fsize/div: 1);
						if (div*nb == fsize)
							G = read_raw_file(s, n, (bigendian->value()? -nb: nb), !bigendian->active());
						}
				}
			}
		}
		if (!G) {
			fl_message("Error loading file!");
			if (vs.size())
				wgl->SetSurface(vs[vs_i]);
			info_box->copy_label("");
		} else {
			if (gridtype) // to set the grid spacing
			{
				accept = 0;
				ratiox->value(1.0); ratioy->value(1.0); ratioz->value(1.0);
				ratiowindow->position(x() + 100, y() + 100);
				ratiowindow->show();
				while(ratiowindow->shown())
					Fl::wait();
				if (accept) {
					G->d[0] = ratiox->value();
					G->d[1] = ratioy->value();
					G->d[2] = ratioz->value();
					for (int i = 0; i != 3; i++)
						G->L[i] = G->d[i]*G->N[i];
				}
			}

			// Init the MC33 temporary structures to use the calc_isosurface function:
			M = create_MC33(G);

			vs.clear();

			numvert_box->copy_label("");
			time_box->copy_label("");
			numvert_box->redraw();
			time_box->redraw();

			// set the scale and center of the surfaces:
			float r[3] = {0.5f*G->L[0], 0.5f*G->L[1], 0.5f*G->L[2]};
#ifndef GRD_orthogonal
			if (G->nonortho)
				mult_Abf(G->_A, r, r, 0);
#endif
			wgl->set_center(r[0] + G->r0[0], r[1] + G->r0[1], r[2] + G->r0[2]);
			int i = G->L[1] > G->L[0];
			if (G->L[2] > G->L[i]) i = 2;
			wgl->set_scale(G->L[i]*0.6f);
			wgl->reset_view();

			// fill the grid info box:
			ostringstream info;
			unsigned int nx = G->N[0] + 1, ny = G->N[1] + 1, nz = G->N[2] + 1;
			info << "** GRID  INFO **\nnum. data:\n" << nx << " x " << ny << " x " << nz;
			info << "\nlengths:\n" << G->L[0] << " x " << G->L[1] << " x " << G->L[2];
			info << "\nratio:\n" << G->d[0] << " : " << G->d[1] << " : " << G->d[2];
			info << "\nmemory used:\n" << (sizeof(_GRD) + nz*(ny + 1)*sizeof(void*) +
			nx*ny*nz*sizeof(GRD_data_type))*9.5367431640625e-7 << " MiB";
			info_box->copy_label(info.str().c_str());
		}
		info_box->redraw();
	}

	static void Delete_cb(Fl_Widget*, void *userdata) {
		MyWindowClass *w = reinterpret_cast<MyWindowClass*>(userdata);
		w->Delete_cb2();
	}

	void Delete_cb2() {
		if (vs.empty())
			return;
		free_surface_memory(vs[vs_i]);
		vs.erase(vs.begin() + vs_i);
		if (vs_i)
			--vs_i;
		fill_surface_info();
	}

	static void Open_cb(Fl_Widget*, void *userdata) {
		MyWindowClass *w = reinterpret_cast<MyWindowClass*>(userdata);
		w->Open_cb2();
	}

	void color_cb2() {
		unsigned char *c = reinterpret_cast<unsigned char*>(&defaultcolor);
		if (fl_color_chooser("Color Chooser", c[0], c[1], c[2])) {
			DefaultColorMC = defaultcolor;
			if (vs.size()) {
				int i = vs[vs_i]->nV;
				while(i--)
					vs[vs_i]->color[i] = defaultcolor;
				wgl->redraw();
			}
		}
	}

	static void color_cb(Fl_Widget*, void *userdata) {
		MyWindowClass *w = reinterpret_cast<MyWindowClass*>(userdata);
		w->color_cb2();
	}

	void save_cb2() {
		if (vs.empty())
			return;
		char *s = fl_file_chooser("Save surface file", "Surface files(*.{txt,sup})", "");
		if (!s)
			return;
		if (ask_replace_file(s)) {
			char *c = strrchr(s, '.');
			if (c) {
				if (strcmp(c, ".sup"))
					write_txt_s(s, vs[vs_i]);
				else
					write_bin_s(s, vs[vs_i]);
			}
		}
	}

	static void save_cb(Fl_Widget*, void *userdata) {
		MyWindowClass *w = reinterpret_cast<MyWindowClass*>(userdata);
		w->save_cb2();
	}

	static void accept_cb(Fl_Widget *w, void *userdata) {
		MyWindowClass *o = reinterpret_cast<MyWindowClass*>(userdata);
		o->accept = 1;
		w->window()->hide();
	}

	static void filedatatype_cb(Fl_Widget *w, void *userdata) {
		MyWindowClass *o = reinterpret_cast<MyWindowClass*>(userdata);
		if (*(w->label()) == 'i')
			o->bigendian->activate();
		else {
			o->bigendian->deactivate();
			o->bigendian->value(0);
		}
		o->reswindow->redraw();
	}

	void ButtonLR_cb(int i) {
		if (vs.empty())
			return;
		if (i) {
			if (++vs_i == vs.size())
				vs_i = 0;
		} else {
			if (vs_i == 0)
				vs_i = vs.size();
			--vs_i;
		}
		fill_surface_info();
	}

	static void ButtonL_cb(Fl_Widget*, void *userdata) {
		MyWindowClass *w = reinterpret_cast<MyWindowClass*>(userdata);
		w->ButtonLR_cb(0);
	}

	static void ButtonR_cb(Fl_Widget*, void *userdata) {
		MyWindowClass *w = reinterpret_cast<MyWindowClass*>(userdata);
		w->ButtonLR_cb(1);
	}

	static void cancel_cb(Fl_Widget *w, void*) {
		w->window()->hide();
	}

public:
	MyWindowClass(int W,int H, const char *L = 0) : Fl_Double_Window(W, H, L) {
		Fl_Group *g = new Fl_Group(W - 110, 0, 110, H - 23);
		isoval_input = new Fl_Value_Input(W - 50, 10, 45, 20, "Isovalue:");
		isoval_input->callback(CalcIso_cb, (void*)this);
		isoval_input->when(FL_WHEN_ENTER_KEY|FL_WHEN_RELEASE);
		Fl_Button *o = new Fl_Button(W - 90, 40, 75, 20,"Calculate");
		o->callback(CalcIso_cb, (void*)this);
		o = new Fl_Button(W - 90, 70, 75, 20,"Color");
		o->callback(color_cb, (void*)this);
		o = new Fl_Button(W - 56, 100, 41, 20,"Save");
		o->callback(save_cb, (void*)this);
		o = new Fl_Button(W - 90, 100, 12, 20,"@<");
		o->callback(ButtonL_cb, (void*)this);
		o = new Fl_Button(W - 78, 100, 12, 20,"@>");
		o->callback(ButtonR_cb, (void*)this);
		o = new Fl_Button(W - 90, 130, 75, 20,"Delete");
		o->callback(Delete_cb, (void*)this);
		o = new Fl_Button(W - 90, 160, 75, 20,"Open grid");
		o->callback(Open_cb, (void*)this);
		info_box = new Fl_Box(W - 106, 190, 106, H - 153);
		info_box->align(FL_ALIGN_TOP|FL_ALIGN_LEFT|FL_ALIGN_INSIDE|FL_ALIGN_WRAP);
		g->resizable(info_box);
		g->end();
		g = new Fl_Group(0, H - 22, W, 22);
		numvert_box = new Fl_Box(0, H - 20, 360, 20);
		numvert_box->box(FL_THIN_DOWN_BOX);
		numvert_box->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
		time_box = new Fl_Box(363, H - 20, W - 363, 20);
		time_box->box(FL_THIN_DOWN_BOX);
		time_box->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
		g->resizable(time_box);
		g->end();
		wgl = new gl_window(0, 2, W - 110, H - 25);
		resizable(wgl);
		end();
		size_range(470, 330);

		reswindow = new Fl_Double_Window(180, 90, "Resolution Nx, Ny, Nz");
			o = new Fl_Button(100, 60, 65, 20, "Accept");
			o->callback(accept_cb, (void*)this);
			res_w[0] = new Fl_Value_Input(11, 10, 46, 20);
			res_w[1] = new Fl_Value_Input(67, 10, 46, 20, ",");
			res_w[2] = new Fl_Value_Input(123, 10, 46, 20, ",");
			Fl_Radio_Round_Button *rb = new Fl_Radio_Round_Button(10, 35, 60, 20, "integer");
			rb->callback(filedatatype_cb, (void*)this);
			rb->value(1);
			rb = new Fl_Radio_Round_Button(10, 60, 60, 20, "real");
			rb->callback(filedatatype_cb, (void*)this);
			bigendian = new Fl_Check_Button(90, 35, 80, 20, "Big endian");
		reswindow->end();
		reswindow->set_modal();

		ratiowindow = new Fl_Double_Window(150, 70, "Aspect ratio X:Y:Z");
			g = new Fl_Group(5, 35, 150, 30, "");
			g->begin();
				o = new Fl_Button(10, 40, 60, 20, "Accept");
				o->callback(accept_cb, (void*)this);
				o = new Fl_Button(80, 40, 60, 20, "Cancel");
				o->callback(cancel_cb);
			g->end();
			ratiox = new Fl_Value_Input(11, 10, 36, 20);
			ratioy = new Fl_Value_Input(57, 10, 36, 20, ":");
			ratioz = new Fl_Value_Input(103, 10, 36, 20, ":");
		ratiowindow->end();
		ratiowindow->set_modal();
		DefaultColorMC = defaultcolor = 0xff669900;
		G = 0;
		M = 0;
	}

	int handle(int event) {
		int ret = Fl_Double_Window::handle(event);
		if (Fl::event_key(FL_Control_L) || Fl::event_key(FL_Control_R)) {
			switch (Fl::event_key()) {
			case FL_Delete: // ctrl + delete, delete the data grid
				free_memory_grd(G);
				G = 0;
				free_MC33(M);
				M = 0;
				info_box->copy_label("");
				info_box->redraw();
				break;
			case 'F':
			case 'f': // ctrl + f, swap back and front faces
				if (vs.size()) {
					{
						int c;
						glGetIntegerv(GL_FRONT_FACE, &c);
						glFrontFace(c == GL_CW? GL_CCW: GL_CW);
					}
					{ // invert the normals
						unsigned int n = 3*vs[vs_i]->nV;
						float *v = vs[vs_i]->N[0];
						for (unsigned int i = 0; i != n; i++)
							v[i] = -v[i];
					}
					wgl->redraw();
				}
				break;
			case 'L':
			case 'l': // ctrl + l, enable/disable back face lighting
				wgl->face_F_FB();
				break;
			}
		}
		return ret;
	}
};

int main(int argc, char **argv) {
	FL_NORMAL_SIZE = 12;
	MyWindowClass mainwin(600, 400, "Marching cubes 33");
	mainwin.show(argc, argv);
	return(Fl::run());
}

