/*
	File: TestMC33_glut.c
	Programmed by: David Vega - dvega@uc.edu.ve
	December 2021
	September 2023
	This is an open source code. The distribution and use rights are under the terms of the MIT license (https://opensource.org/licenses/MIT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef __APPLE__
/* Defined before OpenGL and GLUT includes to avoid deprecation messages */
#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#if defined(_WIN32) && !defined(__CYGWIN__)
#include <io.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif

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

int WIDTH = 600, HEIGHT = 450;

surface** vs = 0;
unsigned int vs_i, vs_n = 0, vs_max = 0;
_GRD *G;
MC33 *MC;

surface *surf;
float bgr, bgg, bgb;
float lightPosition[2][4];
float MGL[16];
float oldM[9], iniV[3];
float tX, tY;
int oldX, oldY;
float scaleSize, scaleFactor;
float cX, cY, cZ;

char surfSizeStr[80], timeStr[32];

char isoString[256], auxInputString[512];
int showMessage, inputIso, saveIso, showHelp, setSurfColor;
float txtF, txtX, txtY;

int drawPoints;
int mouseButton;
int mouseButtonRotate, mouseButtonMove, mouseButtonScale;

void map_to_sphere(float x, float y, float *NewV) {
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

void reset_view() {
	float *f;
	for (f = MGL + 12; --f != MGL;) // MGL[0] is assigned later
		*f = 0.0f;
	for (f = oldM + 8; --f != oldM;) // oldM[0] and oldM[8] are assigned later
		*f = 0.0f;
	oldM[0] = oldM[4] = oldM[8] = MGL[0] = MGL[5] = MGL[10] = scaleFactor = 1.0f;
	tX = tY = 0.0f;
	MGL[15] = scaleSize;
	drawPoints = 0;
	glutPostRedisplay();
}

void reshape(int w, int h) {
	if (w < 120 || h < 90) {
		if (w < 120)
			w = 120;
		if (h < 90)
			h = 90;
		glutReshapeWindow(w, h);
	}
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (w < h) {
		float f = (float)h/w;
		txtF = 8.0f/w;
		txtX = 2*txtF - 1.0f;
		txtY = f - txtF*(26.0f/8);
		glOrtho(-1.0f, 1.0f, -f, f, -10.0f, 10.0f);
	} else {
		float f = (float)w/h;
		txtF = 8.0f/h;
		txtX = 2*txtF - f;
		txtY = 1.0f - txtF*(26.0f/8);
		glOrtho(-f, f, -1.0f, 1.0f, -10.0f, 10.0f);
	}
	glMatrixMode(GL_MODELVIEW);
	WIDTH = w;
	HEIGHT = h;
	glutPostRedisplay(); // update the drawing
}

void set_glMaterial(GLenum face) {
	float mat_specular[4] = { 1.0f, 0.7f, 0.2f, 1.0f};
	float shininess = 50.0f;

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, face == GL_FRONT_AND_BACK);
	glColorMaterial(face, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(face, GL_SPECULAR, mat_specular);
	glMaterialf(face, GL_SHININESS, shininess);
}

void init_gl() {
	float light0Ambient[4] = {0.2f, 0.2f, 0.2f, 1.0f};
	float light0Diffuse[4] = {0.8f, 0.8f, 0.8f, 1.0f};

	float light1Ambient[4] = {0.1f, 0.1f, 0.1f, 1.0f};
	float light1Diffuse[4] = {0.6f, 0.6f, 0.6f, 1.0f};

	// light 0 settings
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0Ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0Diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0Diffuse);

	// light 1 settings
	glEnable(GL_LIGHT1);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light1Ambient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1Diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light1Diffuse);

	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	set_glMaterial(GL_FRONT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(bgr, bgg, bgb, 0.0);
}

void init_param() {
	mouseButtonRotate = GLUT_LEFT_BUTTON;
	mouseButtonMove = GLUT_RIGHT_BUTTON;
	mouseButtonScale = GLUT_MIDDLE_BUTTON;
	tX = tY = bgr = bgg = bgb = cX = cY = cZ = 0.0f;
	scaleFactor = scaleSize = 1.0f;
	lightPosition[0][0] = 1.5f;
	lightPosition[0][1] = -0.5f;
	lightPosition[0][2] = 2.0f;
	lightPosition[0][3] = 0.0f;
	lightPosition[1][0] = -1.5f;
	lightPosition[1][1] = 0.5f;
	lightPosition[1][2] = 2.0f;
	lightPosition[1][3] = 0.0f;
	setSurfColor = inputIso = saveIso = 0;
	showHelp = 1;
	init_gl();
}

// display a char string in the glut window
void displaytext(const char *s, float x, float y) {
	if (isfinite(x))
		glRasterPos2f(x, y);
	while (*s)
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *(s++));
}

void info_text() {
	glLoadIdentity();
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	if (showMessage || saveIso) {
		glColor4f(0.0f, 0.0f, 1.0f, 0.8f);
		float halfw = (saveIso || setSurfColor? 48.0f: 32.0f)*txtF;
		glEnable(GL_BLEND);
		glRectf(-halfw, -16.0f*txtF, halfw, 16.0f*txtF);
		glDisable(GL_BLEND);
		glColor3f(0.75f, 0.75f, 0.75f);
		displaytext("Press Esc to close this box", txtF - halfw, -15.0f*txtF);
		if (saveIso) {
			if (showMessage) {
				displaytext("Type the surface filename (*.txt or *.sup):", -43*txtF, (20.0f/8)*txtF);
				glColor3f(1.0f, 1.0f, 1.0f);
				displaytext(auxInputString, -40*txtF, (-20.0f/8)*txtF);
				glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '_');
			} else
				displaytext("File already exists, overwrite? (y/N)", -37*txtF, 0.0f);
		} else if (setSurfColor) {
			displaytext("Input r (0-255) g (0-255) b (0-255) colors:", -43*txtF, (20.0f/8)*txtF);
			glColor3f(1.0f, 1.0f, 1.0f);
			displaytext(auxInputString, -11*txtF, (-20.0f/8)*txtF);
			glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '_');
		} else
			displaytext(auxInputString, -txtF*strlen(auxInputString), 0.0f);
	} else if (WIDTH > 328 && HEIGHT > 225) {
		if (showHelp) {
			glColor4f(0.20f, 0.20f, 0.20f, 0.6f);
			glEnable(GL_BLEND);
			glRectf(-45.0f*txtF, -27.5f*txtF, 45.0f*txtF, 28.5f*txtF);
			glDisable(GL_BLEND);
			glColor3f(1.0f, 1.0f, 1.0f);
			float f = -41*txtF;
			displaytext("KEYBOARD MENU:", -30*txtF, (90.0f/4)*txtF);
			displaytext("I            : Enter an isovalue", f, (70.0f/4)*txtF);
			displaytext("H            : Show/hide this information", f, (55.0f/4)*txtF);
			displaytext("S            : Save the current surface", f, (40.0f/4)*txtF);
			displaytext("C            : Sets the surface color", f, (25.0f/4)*txtF);
			displaytext("R            : Resets the view", f, (10.0f/4)*txtF);
			displaytext("F            : Swap back and front faces", f, -(5.0f/4)*txtF);
			displaytext("L            : back face lighting on/off", f, -(20.0f/4)*txtF);
			displaytext("Escape key   : Hide help and messages", f, (-35.0f/4)*txtF);
			displaytext("Space bar    : Show the next surface", f, (-50.0f/4)*txtF);
			displaytext("Alt + Space  : Show the previous surface", f, (-65.0f/4)*txtF);
			displaytext("Alt + Delete : Delete the current surface", f, (-80.0f/4)*txtF);
			displaytext("Alt + Q      : Exits the program", f, (-95.0f/4)*txtF);
		} else {
			glColor3f(0.75f, 0.75f, 0.75f);
			displaytext("Press h for help", -32*txtF - txtX , txtY);
		}
	}
		
	glColor3f(0.75f, 0.75f, 0.75f);
	displaytext("Isovalue: ", txtX, txtY);
	if (inputIso) {
		glColor3f(1.0f, 1.0f, 1.0f);
		displaytext(auxInputString, txtX + 20*txtF, txtY);
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '_');
		glColor3f(0.75f, 0.75f, 0.75f);
	} else
		displaytext((surf? isoString: ""), txtX + 20*txtF, txtY);
	displaytext(surfSizeStr, txtX, (-14.0/8)*txtF - txtY);
	//puts(surfSizeStr);
	if (8*(int)(strlen(timeStr) + strlen(surfSizeStr)) + 16 <= WIDTH)
		displaytext(timeStr, -txtX - 2*txtF*strlen(timeStr), (-14.0/8)*txtF - txtY);
	
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
}

// render function
void draw() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition[0]);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPosition[1]);
	MGL[12] = tX*scaleSize - cX*MGL[0] - cY*MGL[4] - cZ*MGL[8];
	MGL[13] = tY*scaleSize - cX*MGL[1] - cY*MGL[5] - cZ*MGL[9];
	MGL[14] = -cX*MGL[2] - cY*MGL[6] - cZ*MGL[10];
	//MGL[15] = scaleFactor*scaleSize;
	glLoadMatrixf(MGL);

	if (surf) {
		if (drawPoints)
			drawdraftsurface(surf);
		else
			drawsurface(surf);
	}

	info_text();
	glutSwapBuffers();
}

// mouse motion
void motion(int x, int y) {
	if (mouseButton == mouseButtonRotate) {
		float qm[9]; // the 3 first elements of qm are the end vector to map to sphere
		// the elements 3 to 6 of qm are the rotation quaternion.
		// To simplify the calculation, the angle of rotation is multiplied by 2.
		map_to_sphere(2.0f*x/(WIDTH - 1) - 1.0f,1.0f - 2.0f*y/(HEIGHT - 1), qm);
		qm[4] = iniV[2]*qm[1] - iniV[1]*qm[2]; // qi
		qm[5] = iniV[0]*qm[2] - iniV[2]*qm[0]; // qj
		qm[6] = iniV[1]*qm[0] - iniV[0]*qm[1]; // qk
		float s = qm[4]*qm[4] + qm[5]*qm[5] + qm[6]*qm[6]; // sin^2(angle)
		if (s < 0.00001f*0.00001f)
			return;
		qm[3] = iniV[0]*qm[0] + iniV[1]*qm[1] + iniV[2]*qm[2]; // qr, cos(angle) instead cos(angle/2)
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
	} else if (mouseButton == mouseButtonMove) {
		float t = 2*(WIDTH > HEIGHT ? scaleFactor/HEIGHT : scaleFactor/WIDTH);
		tX += (x - oldX)*t;
		tY -= (y - oldY)*t;
		oldX = x;
		oldY = y;
	} else if (mouseButton == mouseButtonScale) {
		scaleFactor += (y - oldY)*scaleFactor*0.001f;
		MGL[15] = scaleFactor*scaleSize;
		oldY = y;
	}
	glutPostRedisplay();
}

void show_message(const char *msg) {
	showMessage = 1;
	strncpy(auxInputString, msg, sizeof auxInputString - 1);
	glutMouseFunc(0);
	glutMotionFunc(0);
	glutPostRedisplay();
}

void mouse(int, int, int, int);
void hide_message() {
	showMessage = 0;
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
}

void fill_surface_info() {
	if (!vs_n) {
		surfSizeStr[0] = 0;
		timeStr[0] = 0;
		surf = 0;
	} else {
		int n;
		surf = vs[vs_i];
		sprintf(surfSizeStr, "Surface: %d|%d  Vertices: %d  Triangles: %d", vs_i + 1, vs_n, surf->nV, surf->nT);
		sprintf(timeStr, "Time: %d ms", surf->user.i[0]);
		sprintf(isoString, "%.6f", surf->user.df[1]);
		n = strlen(isoString);
		while (--n > 0) {
			if (isoString[n] == '0')
				isoString[n] = (char)0;
			else
				break;
		}
		if (isoString[n] == '.')
			isoString[n] = (char)0;
	}
	glutPostRedisplay();
}

int file_exist(const char *filename) {
	return (access(filename, F_OK)? 0: 1);
}

void delete_iso() {
	unsigned int i;
	if (!vs_n)
		return;
	free_surface_memory(vs[vs_i]);
	--vs_n;
	for (i = vs_i; i < vs_n; i++)
		vs[i] = vs[i + 1];
	if (vs_i)
		--vs_i;
	fill_surface_info();
}

void calc_iso(double iso) {
	unsigned int i;
	for (i = 0; i != vs_n; i++)
		if (vs[i]->user.df[1] == iso) {
			vs_i = i;
			fill_surface_info();
			return;
		}
	if (vs_max == vs_n) {
		vs_max += 16;
		void *t = realloc(vs, vs_max*sizeof(void*));
		if (t)
			vs = (surface**)t;
		else {
			vs_max = vs_n;
			show_message("Insufficient memory!");
			return;
		}
	}
	if (MC) {
		int t = clock();
		surface *S = calculate_isosurface(MC, iso);
		t = clock() - t;
		if (S) {
			if (S->nV > 0) {
				S->user.i[0] = t*1000/CLOCKS_PER_SEC;
				S->user.df[1] = iso; // save the isovalue as double for better comparison.
				vs_i = vs_n++;
				vs[vs_i] = S;
				fill_surface_info();
			} else {
				show_message("Empty surface!");
				free_surface_memory(S);
			}
		} else
			show_message("Insufficient memory!");
	} else
		show_message("No grid data!");
}

void select_iso(int next) {
	if (!vs_n)
		return;
	if (next) {
		if (++vs_i == vs_n)
			vs_i = 0;
	} else {
		if (vs_i == 0)
			vs_i = vs_n;
		--vs_i;
	}
	fill_surface_info();
}

void save_iso2() {
	const char *c = strrchr(auxInputString, '.');
	if (c) {
		if (strcmp(c,".sup"))
			c = 0;
	}
	if (c)
		write_bin_s(auxInputString, surf);
	else
		write_txt_s(auxInputString, surf);
}

// called when a mouse button is pressed
void mouse(int button, int state, int x, int y) {
	if (state == GLUT_DOWN) { // if any button is pressed
		oldX = x;
		oldY = y;
		if (button == mouseButtonRotate)
			map_to_sphere(2.0f*oldX/(WIDTH - 1) - 1.0f,1.0f - 2.0f*oldY/(HEIGHT - 1), iniV);
		mouseButton = button;
		drawPoints = 1;
	} else {
		memcpy(oldM, MGL, 3*sizeof(float));
		memcpy(oldM + 3, MGL + 4, 3*sizeof(float));
		memcpy(oldM + 6, MGL + 8, 3*sizeof(float));
		drawPoints = 0;
		glutPostRedisplay();
	}
}

void keyboard(unsigned char c, int, int);

void keyb_string(unsigned char c, int x, int y) {
	int n = strlen(auxInputString);
	if (c >= ' ' && c <= '}' && n + 1 < (int)(sizeof auxInputString)) {
		auxInputString[n] = (char)c;
		auxInputString[n + 1] = (char)0;
	} else if (c == '\b' && n)
		auxInputString[n - 1] = (char)0;
	else if (c == '\r' || c == 27) {
		if (inputIso) {
			if (c == '\r') {
				auxInputString[sizeof isoString - 1] = (char)0;
				double iso = strtod(auxInputString, 0);
				strcpy(isoString, auxInputString);
				calc_iso(iso);
			}
			inputIso = 0;
		} else if (saveIso) {
			hide_message();
			if (auxInputString[0] == (char)0)
				saveIso = 0;
			else if (!file_exist(auxInputString)) {
				saveIso = 0;
				save_iso2();
			}
		} else if (setSurfColor) {
			setSurfColor = 0;
			hide_message();
			if (n > 4) {
				int rgb[3];
				union {
					unsigned char c[4];
					unsigned int ui;
				} u;
				u.ui = 0xff000000;
				sscanf(auxInputString, "%d %d %d", rgb, rgb + 1, rgb + 2);
				for (int i = 0; i != 3; i++)
					u.c[i] = rgb[i]&0xff;
				DefaultColorMC = u.ui;
				if (surf) {
					int i = surf->nV;
					while(i--)
						surf->color[i] = u.ui;
				}
			}
		}
		glutKeyboardFunc(keyboard);
	}
	glutPostRedisplay();
}

void iso_input() {
	auxInputString[0] = (char)0;
	glutKeyboardFunc(keyb_string);
	inputIso = 1;
	hide_message();
}

void set_color() {
	auxInputString[0] = (char)0;
	glutKeyboardFunc(keyb_string);
	showMessage = setSurfColor = 1;
}

void save_iso() {
	if (!vs_n)
		return;
	auxInputString[0] = (char)0;
	saveIso = showMessage = 1;
	glutKeyboardFunc(keyb_string);
}

void keyboard(unsigned char key, int x, int y) {
	if (saveIso) {
		if (key == 'Y' || key == 'y')
			save_iso2();
		saveIso = 0;
	} else
		switch(key) {
		case 27: // escape key
			showHelp = 0;
			hide_message();
			break;
		case ' ':
			select_iso(GLUT_ACTIVE_ALT != glutGetModifiers());
			break;
		case '\b': // Backspace key
		case 127: // Delete key
			if (GLUT_ACTIVE_ALT == glutGetModifiers()) // if the Control key is pressed
				delete_iso();
			break;
		case 'C':
		case 'c':
			set_color();
		break;
		case 'F':
		case 'f': // alt + f, swap back and front faces
			if (vs_n) {
				{
					int c;
					glGetIntegerv(GL_FRONT_FACE, &c);
					glFrontFace(c == GL_CW? GL_CCW: GL_CW);
				}
				{ // invert the normals
					unsigned int n = 3*surf->nV;
					float *v = surf->N[0];
					for (unsigned int i = 0; i != n; i++)
						v[i] = -v[i];
				}
			}
			break;
		case 'H':
		case 'h':
			if (showMessage) {
				showHelp = 1;
				hide_message();
			} else
				showHelp = !showHelp;
			break;
		case 'I':
		case 'i':
			iso_input();
			break;
		case 'L':
		case 'l': // alt + l, enable/disable back face lighting
			{
				GLboolean v;
				glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE, &v);
				set_glMaterial(v? GL_FRONT: GL_FRONT_AND_BACK);
			}
			break;
		case 'Q':
		case 'q':
			if (GLUT_ACTIVE_ALT == glutGetModifiers()) // if the Control key is pressed
				exit(0);
			break;
		case 'R':
		case 'r':
			reset_view();
			break;
		case 'S':
		case 's':
			save_iso();
			break;
		}
	glutPostRedisplay();
}

void terminal_input(char *s, const char *msg) {
	puts(msg);
	fgets(s, 255, stdin);
}

int input_res(int *n, long long size) {
	if (size < 8)
		return 0;
	char s[256];
	terminal_input(s, "Grid dimensions (Nx Ny Nz): ");
	if (strlen(s) > 4) {
		sscanf(s, "%d %d %d", n, n + 1, n + 2);
	} else
		return 0;
	long long div = n[0]*n[1]*n[2];
	if (!div)
		return 0;
	n[5] = n[3] = size/div;
	if (n[3] == 4)
		terminal_input(s, "Data type r (real) or i (integer): ");
	else {
		s[0] = (n[3] < 4? 'i': 'r');
	}
	switch (s[0]) {
	case 'R':
	case 'r':
		n[4] = 1;
		if (n[3] != 4 && n[3] != 8)
			return 0;
		break;
	case 'I':
	case 'i':
		n[4] = 0;
		if (n[3] != 1 && n[3] != 2 && n[3] != 4)
			return 0;
		if (n[3] > 1) {
			terminal_input(s, "big-endian (b) or little-endian (l): ");
			if (s[0] == 'b' || s[0] == 'B')
				n[5] = -n[5];
		}
		break;
	default:
		return 0;
	}
	return div*n[3] == size;
}

long long filesize(const char *filename) {
	long long size = 0;
	FILE *f = fopen(filename, "rb");
	if (f) {
		fseek(f, 0, SEEK_END);
		size = ftell(f);
		fclose(f);
	}
	return size;
}

void set_grid_spacing(double *r, int setr) {
	if (setr < 0) {
		char s[256];
		terminal_input(s, "Grid spacings (dx dy dz): ");
		if (strlen(s) > 4) {
			sscanf(s, "%lf %lf %lf", r, r + 1, r + 2);
			for (int i = 0; i != 3; i++)
				if (r[i] <= 0.0)
					r[i] = 1.0;
		}
	}
}


// output the grid info in console terminal:
void display_grid_info() {
	unsigned int nx = G->N[0] + 1, ny = G->N[1] + 1, nz = G->N[2] + 1;
	printf("\n ***** GRID INFO *****\nNum. data:\n\t%d x %d x %d", nx, ny, nz);
	printf("\nDimensions:\n\t%f x %f x %f", G->L[0], G->L[1], G->L[2]);
	printf("\nOrigin:\n\t(%f, %f, %f)", G->r0[0], G->r0[1], G->r0[2]);
	printf("\nstep:\n\t%f, %f, %f", G->d[0], G->d[1], G->d[2]);
	printf("\nmemory used:\n\t%f MiB\n\n", (sizeof(_GRD) + nz*(ny + 1)*sizeof(void*) +
	nx*ny*nz*sizeof(GRD_data_type))*9.5367431640625e-7);
}


int open_grid_file(const char *s) {
	int setr = -1;
	double r[3] = {1.0, 1.0, 1.0};
	const char *c = strrchr(s, '.');
	if (c) {
		if (!strcmp(c,".grb"))
			G = read_grd_binary(s);
		else if (!strcmp(c,".grd"))
			G = read_grd(s);
		else if (!strcmp(c,".dat")) {
			set_grid_spacing(r, setr);
			setr = 1;
			G = read_dat_file(s);
		}
	}
	if (!G) {
		long long size = filesize(s);
		switch ((int)size) {
		case 8192:
			set_grid_spacing(r, setr);
			G = read_scanfiles(s, 64, 0);
			break;
		case 131072:
			set_grid_spacing(r, setr);
			G = read_scanfiles(s, 256, 1);
			break;
		case 524288:
			set_grid_spacing(r, setr);
			G = read_scanfiles(s, 512, 1);
			break;
		default:
			{
				int n[6];
				if (input_res(n, size)) {
					set_grid_spacing(r, setr);
					G = read_raw_file(s, (unsigned int*)n, n[5], n[4]);
				}
			}
		}
	}
	if (!G) {
		puts("\nError loading grid file!\n");
	} else {
		surf = 0;
		if (setr)
			for (int i = 0; i != 3; i++) {
				G->d[i] = r[i];
				G->L[i] = r[i]*G->N[i];
			}
	}
	return G != 0;;
}

// set the scale and center of the graphic window
void set_scale() {
	float r[3] = {0.5f*G->L[0], 0.5f*G->L[1], 0.5f*G->L[2]};
#ifndef GRD_orthogonal
	if (G->nonortho)
		mult_Abf(G->_A, r, r, 0);
#endif
	cX = r[0] + G->r0[0];
	cY = r[1] + G->r0[1];
	cZ = r[2] + G->r0[2];
	int i = G->L[1] > G->L[0];
	if (G->L[2] > G->L[i])
		i = 2;
	scaleSize = (2.0f/3)*G->L[i];
}

/*********
Test functions used to generate the grid (implicit surfaces)
double f(double x, double y, double z)
**********/

double f1(double x, double y, double z) {
	return (sin(x*y) + sin(y*z) + sin(x*z))/(1.0 + x*x + y*y + z*z);
}

// https://blogs.ams.org/visualinsight/2016/04/15/barth-sextic/
#define phi ((sqrt(5) + 1)/2)
double Barthsextic(double x, double y, double z) {
	const double phi2 = phi*phi, w2 = 1.0*1.0;
	x *= x; y *= y; z *= z;
	double t = (x + y + z - w2)*w2;
	return 4*(phi2*x - y)*(phi2*y - z)*(phi2*z - x) - (1 + 2*phi)*t*t;
}

double ga(double cx, double cy, double cz, double x, double y, double z) {
	x -= cx; y -= cy; z -= cz;
	return exp(-0.5*(x*x + y*y + z*z));
}

double fourGaussian(double x, double y, double z) {
	return ga(-1.5, -1.5, -1.5, x, y, z) + ga(1.5, 1.5, -1.5, x, y, z) +
		ga(1.5, -1.5, 1.5, x, y, z) + ga(-1.5, 1.5, 1.5, x, y, z);
}

double Cube_f(double x, double y, double z) {
	x = fabs(x); y = fabs(y); z = fabs(z);
	if (x > y) {
		if (x > z)
			z = x;
	} else {
		if (y > z)
			z = y;
	}
	return 1.0 - (1.0/3.0)*z;
}

double Cylinder_f(double x, double y, double z) {
	double t = sqrt(x*x + z*z);
	y = sqrt(2.0)*fabs(y);
	if (t < y)
		t = y;
	return 1.0 - (1.0/3.0)*t;
}

//(9z^2 - 1)(1 - z^2) - 2y(y^2 - 3x^2)(1 - z^2) - (x^2+y^2)^2
// https://en.wikipedia.org/wiki/Implicit_surface
double Genus2(double x, double y, double z) {
	double t = y*y;
	x *= x; z = z*z - 1.0;
	y = (2.0*y*(t - 3.0*x) - 9.0*z - 8.0)*z;
	x += t;
	return y - x*x;
}

// -x^4 + 5*x^2 - y^4 + 5*y^2 - z^4 + 5*z^2 - 11.8
// https://www-sop.inria.fr/galaad/surface/
double TangleCube(double x, double y, double z) {
	x *= x; y *= y; z *= z;
	return x*(5 - x) + y*(5 - y) + z*(5 - z) - 11.8;
}

// https://en.wikipedia.org/wiki/Implicit_surface#Smooth_approximations_of_several_implicit_surfaces
double threeTori(double x, double y, double z) {
	const double R2 = 1.0*1.0, a2 = 0.2*0.2, r = 0.01;
	x *= x; y *= y; z *= z;
	double t = x + y + z + R2 - a2;
	t *= t;
	return r - (t - 4*R2*(x + y))*(t - 4*R2*(x + z))*(t - 4*R2*(y + z));
}

// http://virtualmathmuseum.org/Surface/decocube/decocube.html
const double cc = 1.3;
double DecoCube(double x, double y, double z) {
	const double c2 = 2.0 - cc*cc, ff = 0.02;
	x = x*x - 1.0; y = y*y - 1.0; z = z*z - 1.0;
	double t1 = x + y + c2, t2 = y + z + c2, t3 = z + x + c2;
	t1 *= t1; t2 *= t2; t3 *= t3;
	return ff - (t1 + z*z)*(t2 + x*x)*(t3 + y*y);
}

// https://math.stackexchange.com/q/46222
double f(double x, double y, double z) {
	x *= x; y *= y; z *= z;
	double t1 = (2.92*(x - 1)*x + 1.7*y)*(y - 0.88);
	double t2 = (2.92*(y - 1)*y + 1.7*z)*(z - 0.88);
	double t3 = (2.92*(z - 1)*z + 1.7*x)*(x - 0.88);
	return 0.04 - t1*t1 - t2*t2 - t3*t3;
}
//(2.92(x-1)x^2(x+1)+1.7y^2)^2*(y^2-0.88)^2 +
//(2.92(y-1)y^2(y+1)+1.7z^2)^2*(z^2-0.88)^2 +
//(2.92(z-1)z^2(z+1)+1.7x^2)^2*(x^2-0.88)^2 - 0.04
/******* End of test functions *********/

int main(int argc, char **argv) {

	int i = 1;
	if (argc == 1 || !open_grid_file(argv[1])) {
		printf("To read a grid file, close this program and type:\n%s filename\n", argv[0]);
		puts("\nSelect one of the default grids:");
		puts("\t1. (sin(x*y) + sin(y*z) + sin(x*z))/(1 + x*x + y*y + z*z)");
		puts("\t2. Barth sextic");
		puts("\t3. Four Gaussian functions, isovalue (0, 1)");
		puts("\t4. Cube, isovalue [0, 1)");
		puts("\t5. Cylinder, isovalue [-0.4142, 1)");
		puts("\t6. (9z^2 - 1)(1 - z^2) - 2y(y^2 - 3x^2)(1 - z^2) - (x^2+y^2)^2");
		puts("\t7. -x^4 + 5*x^2 - y^4 + 5*y^2 - z^4 + 5*z^2 - 11.8");
		puts("\t8. Three intersecting tori, isovalue 0");
		puts("\t9. DecoCube, isovalue 0");
		puts("\t0. Leo cube?, isovalue 0");
		do
			i = getchar();
		while (!i);
		switch (i) {
			case '1':
				G = generate_grid_from_fn(-3.0, -3.0, -3.0, 3.0, 3.0, 3.0, 0.04, 0.04, 0.04, f1);
				break;
			case '2':
				G = generate_grid_from_fn(-2.0, -2.0, -2.0, 2.0, 2.0, 2.0, 0.02, 0.02, 0.02, Barthsextic);
				break;
			case '3':
				G = generate_grid_from_fn(-4.0, -4.0, -4.0, 4.0, 4.0, 4.0, 0.04, 0.04, 0.04, fourGaussian);
				break;
			case '4':
				G = generate_grid_from_fn(-3.0, -3.0, -3.0, 3.0, 3.0, 3.0, 0.04, 0.04, 0.04, Cube_f);
				break;
			case '5':
				G = generate_grid_from_fn(-3.0, -3.0, -3.0, 3.0, 3.0, 3.0, 0.04, 0.04, 0.04, Cylinder_f);
				break;
			case '6':
				G = generate_grid_from_fn(-1.875, -2.2, -1.875, 1.875, 1.55, 1.875, 0.025, 0.025, 0.025, Genus2);
				break;
			case '7':
				G = generate_grid_from_fn(-3.0, -3.0, -3.0, 3.0, 3.0, 3.0, 0.03, 0.03, 0.03, TangleCube);
				break;
			case '8':
				G = generate_grid_from_fn(-1.5, -1.5, -1.5, 1.5, 1.5, 1.5, 0.015, 0.015, 0.015, threeTori);
				break;
			case '9':
				G = generate_grid_from_fn(-1.5, -1.5, -1.5, 1.5, 1.5, 1.5, 0.015, 0.015, 0.015, DecoCube);
				break;
			case '0':
				G = generate_grid_from_fn(-1.2, -1.2, -1.2, 1.2, 1.2, 1.2, 0.006, 0.006, 0.006, f);
				break;
			default:
				i = 0;
		}
	}
	if (i) {
		display_grid_info();
		i = 1;
		glutInit(&i, argv);
		glutInitWindowSize(WIDTH, HEIGHT);
		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
		glutInitWindowPosition(50, 50);
		glutCreateWindow("Marching cubes 33");

		init_param();
		MC = create_MC33(G);
		set_scale();
		reset_view();

		glutDisplayFunc(draw); // set the render function
		glutKeyboardFunc(keyboard); // register the keyboard functions
		glutReshapeFunc(reshape); // register the reshape function
		hide_message();
		glutMainLoop();
	}
	return 0;
}

