#include "vis.h"

#define THROWERROR(msg) { printf("ERROR: %s(): %s line %d: \"%s\"\n", __func__, __FILE__, __LINE__, msg); exit(-1); }


vis::vis(int _target_fps, int argc, char **argv) : 
	GLWindow((char *)"OpenGL", 10, 10, 1024, 768, GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH, argc, argv) {
	int i, j;

	factor = 0.120;

	framecount = 0;
	target_fps = _target_fps;
	vblank_sync = false;
	data_invalidated = true;
	view_invalidated = true;
	
	r = NULL;
	rad = NULL;
	types = NULL;

	/* move and display timers as close to target fps as possible */
	display_timer_id = SetTimer(1000/target_fps, true);
	move_timer_id = SetTimer(1000/target_fps, true);
	/* info every 2 seconds */
	info_timer_id = SetTimer(1000, true);
	printf( "Simple OpenGL visualisation window; press '?' to display help and 'q' to quit.\n" );
	
	base_cam->location.z = -200.0;

	base_cam->type = JCAMERA_TYPE_WATCHORIGIN;
	

	base_light->pos.z = -1000.0; //-5.0;
//	base_light->zrot(M_PI/4.0);
	base_light->yrot(-M_PI/6.0);

	GLfloat colors[3] = { 1.0f, 1.0f, 1.0f };
	glClearColor(colors[0], colors[1], colors[2], 1.0f);

// Fog code
//	GLfloat fog_color[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
//	glFogi(GL_FOG_MODE, GL_LINEAR);
//	glFogfv(GL_FOG_COLOR, fog_color);
//	glHint(GL_FOG_HINT, GL_DONT_CARE); 
//	glFogf(GL_FOG_START, 160.0); 
//	glFogf(GL_FOG_END, 400.0);
//	glFogf(GL_FOG_DENSITY, 1.0f);
//	glEnable(GL_FOG);


//			glDepthFunc(GL_LEQUAL); 
//			glEnable(GL_DEPTH_TEST); 
//			glShadeModel(GL_FLAT);

	glutPostRedisplay();

// Set up colors for the system beads
	ran1_seed = -1;
	for (i = 0; i < MAX_TYPES; i++) {

		switch (i) {
			case 0: {
				bead_amb[(i*4)] = 0.6f;
				bead_amb[(i*4)+1] = 0.6f;
				bead_amb[(i*4)+2] = 0.6f;
				bead_amb[(i*4)+3] = 1.0f;
				break;
			}
			case 1: {
				bead_amb[(i*4)] = 0.2f;
				bead_amb[(i*4)+1] = 0.5f;
				bead_amb[(i*4)+2] = 0.5f;
				bead_amb[(i*4)+3] = 1.0f;
				break;
			}
			case 2: {
				bead_amb[(i*4)] = 0.2f;
				bead_amb[(i*4)+1] = 0.2f;
				bead_amb[(i*4)+2] = 0.5f;
				bead_amb[(i*4)+3] = 1.0f;
				break;
			}
			case 3: {
				bead_amb[(i*4)] = 0.5f;
				bead_amb[(i*4)+1] = 0.2f;
				bead_amb[(i*4)+2] = 0.2f;
				bead_amb[(i*4)+3] = 1.0f;
				break;
			}
			case 4: {
				bead_amb[(i*4)] = 0.2f;
				bead_amb[(i*4)+1] = 0.5f;
				bead_amb[(i*4)+2] = 0.2f;
				bead_amb[(i*4)+3] = 1.0f;
				break;
			}
			case 5: {
				bead_amb[(i*4)] = 0.4f;
				bead_amb[(i*4)+1] = 0.1f;
				bead_amb[(i*4)+2] = 0.4f;
				bead_amb[(i*4)+3] = 1.0f;
				break;
			}
			case 6: {
				bead_amb[(i*4)] = 0.4f;
				bead_amb[(i*4)+1] = 0.4f;
				bead_amb[(i*4)+2] = 0.1f;
				bead_amb[(i*4)+3] = 1.0f;
				break;				
			}
			default: {
				bead_amb[(i*4)] = 0.7f;
				bead_amb[(i*4)+1] = 0.0f;
				bead_amb[(i*4)+2] = 0.7f;
				bead_amb[(i*4)+3] = 1.0f;
				break;				
			}
		}

		
		for(j=0; j < 4; j++)
		{
			if (j < 3) {	
				bead_diff[(i*4)+j] = bead_amb[(i*4)+j]*0.9;
				bead_spec[(i*4)+j] = bead_amb[(i*4)+j]*0.3f;
				bead_emit[(i*4)+j] = bead_amb[(i*4)+j]*0.3f;
			} else {
				bead_diff[(i*4)+j] = bead_amb[(i*4)+j];
				bead_spec[(i*4)+j] = bead_amb[(i*4)+j];
				bead_emit[(i*4)+j] = bead_amb[(i*4)+j];
			}
		}
		bead_shiny[i] = 25.0;
		show_types[i] = true;
	}
	
	cell_col[0] = 100.0; cell_col[1] = 100.0; cell_col[2] = 100.0; // luke
		
	scene_dlist = glGenLists(1);
	glNewList( scene_dlist, GL_COMPILE );
	glEndList();

// display list for a bead
//	bead_dlist = glGenLists( 1 );
//	glNewList( bead_dlist, GL_COMPILE );
//		glutSolidSphere( 2.0, 50, 50 );
//	glEndList();
}

vis::~vis()
{
	if (types != NULL) free(types);
	if (r != NULL) free(r);
	if (rad != NULL) free(rad);
	glDeleteLists(scene_dlist, 1);
//	glDeleteLists(bead_dlist, 1);
}


void vis::Timer(GLint value) {
	
	float rotspeed = 0.020;
	
	if (value == display_timer_id ) Display();
	else if (value == move_timer_id) {
		if (special_keys[GLUT_KEY_UP] == 1) { 
			if (ORTHOGRAPHIC) {
				if (factor <= 0.004) factor -= 0.001;
				else factor -= 0.004;
				if (factor <= 0.0) factor = 0.001;
			} else for (int i = 0; i < 5; i++) base_cam->location += base_cam->look;
			view_invalidated = true; 
//			_glutReshape(windowdims[2], windowdims[3]);
		}
		else if (special_keys[GLUT_KEY_DOWN] == 1) {
			if (ORTHOGRAPHIC) {
				if (factor < 0.004) factor += 0.001;
				else factor += 0.004;
			} else for (int i = 0; i < 5; i++) base_cam->location -= base_cam->look;
			view_invalidated = true;
//			glOrtho(-windowdims[2]*factor, windowdims[2]*factor, -windowdims[3]*factor, windowdims[3]*factor, 0.01, 1000000.0);
//			_glutReshape(windowdims[2], windowdims[3]);
		}

		if (special_keys[GLUT_KEY_LEFT] == 1) {
			base_cam->zrot(rotspeed); base_light->zrot(rotspeed); 
			view_invalidated = true; 
		}
		else if (special_keys[GLUT_KEY_RIGHT] == 1) { 
			base_cam->zrot(-rotspeed); base_light->zrot(-rotspeed); 
			view_invalidated = true; 
		}
	}
	else if (value == info_timer_id) {
		sprintf(frame_string, "%.1f fps", (1000.0/timers[info_timer_id][0])*(framecount+1));
		view_invalidated = 1;
		Display();
		framecount = 0;
	}
}


void ViewOrtho(int width, int height, float factor) {
	glMatrixMode(GL_PROJECTION);			// Select Projection
	glLoadIdentity();						// Reset The Matrix
	glOrtho(-width*factor, width*factor, -height*factor, height*factor, 0.01, 1000000.0);
	glMatrixMode(GL_MODELVIEW);				// Select Modelview Matrix
	glLoadIdentity();						// Reset The Matrix
}



void vis::Display()
{
	int i, j;
	
	if (ShouldQuitThread() || (!view_invalidated && !data_invalidated)) return;

	if (ORTHOGRAPHIC) ViewOrtho(windowdims[2], windowdims[3], factor);
			
	glLoadIdentity();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	base_cam->see();
	base_light->shine();


	if (view_invalidated) view_invalidated = false;

	if (data_invalidated) {


		
		glDeleteLists(scene_dlist, 1);
		scene_dlist = glGenLists(1);
		glNewList(scene_dlist, GL_COMPILE);

		// draw cell; switch off lighting or funky effects will ensue from cell vertices.

		double x, y, z;

		glDisable(GL_LIGHTING);
			x = cell[0];
			y = cell[1];
			z = cell[2];
			// cell color
			glColor3f(0.0, 0.0, 0.0);
			glEnable (GL_LINE_SMOOTH);
			glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
			glLineWidth(3.0);
			// front cell wall
			glBegin(GL_LINE_STRIP);
				glVertex3f(-x/2.0, -y/2.0, -z/2.0);
				glVertex3f(-x/2.0,  y/2.0, -z/2.0);
				glVertex3f(x/2.0,  y/2.0, -z/2.0);
				glVertex3f(x/2.0, -y/2.0, -z/2.0);
				glVertex3f(-x/2.0, -y/2.0, -z/2.0);
			glEnd();
				
			// rear cell wall
			glBegin(GL_LINE_STRIP);
				glVertex3f(-x/2.0, -y/2.0, z/2.0);
				glVertex3f(-x/2.0,  y/2.0, z/2.0);
				glVertex3f(x/2.0,  y/2.0, z/2.0);
				glVertex3f(x/2.0, -y/2.0, z/2.0);
				glVertex3f(-x/2.0, -y/2.0, z/2.0);
			glEnd( );	
			// right cell wall
			glBegin(GL_LINE_STRIP);
				glVertex3f(-x/2.0, -y/2.0, -z/2.0);
				glVertex3f(-x/2.0,  y/2.0, -z/2.0);
				glVertex3f(-x/2.0,  y/2.0,  z/2.0);
				glVertex3f(-x/2.0, -y/2.0,  z/2.0);
				glVertex3f(-x/2.0, -y/2.0, -z/2.0);
			glEnd();	
			// left cell wall
			glBegin(GL_LINE_STRIP);
				glVertex3f(x/2.0, -y/2.0, -z/2.0);
				glVertex3f(x/2.0,  y/2.0, -z/2.0);
				glVertex3f(x/2.0,  y/2.0,  z/2.0);
				glVertex3f(x/2.0, -y/2.0,  z/2.0);
				glVertex3f(x/2.0, -y/2.0, -z/2.0);
			glEnd();	
		glEnable(GL_LIGHTING);
//		glDisable(GL_BLEND);

		/* Sites */
		j = -1;
		if( nsites > 0 && r != NULL && types != NULL && rad != NULL )
		{
			for( i=0; i<nsites; i++ )
			{
				if( show_types[ types[i] ] == true )
				{
					/*
						Do we need to switch material? Should really sort sites by type, and hence minimise context swiches!
					*/
					if( types[i] != j )
					{
						j = types[i];
						glMaterialf(GL_FRONT, GL_SHININESS, bead_shiny[j]);
						glMaterialfv(GL_FRONT, GL_AMBIENT, &bead_amb[j*4]);
						glMaterialfv(GL_FRONT, GL_DIFFUSE, &bead_diff[j*4]);
						glMaterialfv(GL_FRONT, GL_SPECULAR, &bead_spec[j*4]);
						glMaterialfv(GL_FRONT, GL_EMISSION, &bead_emit[j*4]);
					}
					glPushMatrix();
					// TAGGED SECTION
					glTranslatef(r[(i*3)+0], r[(i*3)+1], r[(i*3)+2]);
					glutSolidSphere(rad[i], 10, 10);
					glPopMatrix();
				}
			}
		}
		

		glEndList();
		data_invalidated = false;
		BaseThreadUnlock();
	}

	glCallList(scene_dlist);

	glFlush();
	glutSwapBuffers();
	
	framecount++;
}


void vis::MouseButton(GLint button, GLint state, GLint x, GLint y) {
	lastx = x;
	lasty = y;
}



void vis::MouseMotion(GLint x, GLint y) {
	if (lastx-x != 0) {
		base_cam->yrot((float)(lastx-x)/windowdims[2] * 2.0);
		base_light->yrot((float)(lastx-x)/windowdims[2] * 2.0);
	}
	if (lasty-y != 0) {
		base_cam->xrot(-(float)(lasty-y)/windowdims[3] * 2.0);
		base_light->xrot(-(float)(lasty-y)/windowdims[3] * 2.0);
	}
	
	lastx = x;
	lasty = y;
	view_invalidated = true;
}



void vis::KeyDown(GLint key, GLint x, GLint y) {
#ifdef __APPLE__
	const GLint vblank_state[] = { 0, 1 };
#endif	

	if (key == 'v') {
		vblank_sync = !vblank_sync;
		#ifdef __APPLE__
			CGLSetParameter(CGLGetCurrentContext(), kCGLCPSwapInterval, &vblank_state[vblank_sync]);
		#endif
		printf("VBlank set to %d\n", vblank_sync);
	}
	else if (key == 's') { base_light->show_pos = !base_light->show_pos; data_invalidated = true; }
	else if (key == 'q') { Stop(); }
	else if (key >= '1' && key <= '9') { show_types[key - '1'] = !show_types[key - '1']; data_invalidated = true; }
}



void vis::Reshape(GLint width, GLint height) {
	view_invalidated = true;
}


void vis::setLightPos(GLfloat x, GLfloat y, GLfloat z) {
	BaseThreadLock();
		base_light->pos.x = x;
		base_light->pos.y = y;
		base_light->pos.z = z;
		data_invalidated = true;
	BaseThreadUnlock();
}

void vis::Update(int _n, double * _r, int * _t, double * _cell, double *_rad) {
	if ((_r == NULL) || (_t == NULL) || (_rad == NULL) || (_n < 0)) return;
	
	BaseThreadLock();
		if (nsites != _n)
		{
			nsites = _n;
			if (r != NULL) free(r);
			if (types != NULL) free(types);
			if (rad != NULL) free(rad);
			if (_n == 0) { 
				r = NULL; 
				types = NULL;
				rad = NULL; 
				return; 
			}
			types = (int *)malloc(sizeof(int)*nsites);
			r = (double *)malloc(sizeof(double)*3*nsites);
			rad = (double *)malloc(sizeof(double)*nsites);
			if(types == NULL || r == NULL || rad == NULL)
			{
				printf("Error: Unable to reallocate types, r, or rad array.\n");
				exit(-1);
			}
		}
		memcpy(r, _r, sizeof(double)*3*nsites);
		memcpy(types, _t, sizeof(int)*nsites);
		memcpy(rad, _rad, sizeof(double)*nsites);
		memcpy(cell, _cell, sizeof(double)*3);
		view_invalidated = true;
		data_invalidated = true;
	BaseThreadUnlock();
}





int port;
	
void *visFunc(void *windowptr)
{
	vis *wnd = (vis *)windowptr;
	jsocket *js;
	
	int connection;
	struct sockaddr_storage caddr;
	socklen_t clen;

	double cell[] = { 10.0, 10.0, 10.0 };
	int n, new_n, *types;
	double *r, *rad;
	
	if (wnd == NULL)
	{
		printf("Error: NULL window pointer paseed to %s\n", __func__);
		exit(-1);
	}
	
	n = 1;
	types = (int *)malloc(sizeof(int)*n);
	r = (double *)malloc(sizeof(double)*n*3);
	rad = (double *)malloc(sizeof(double)*n);
	
	types[0] = 0;
	r[0] = 0.0;
	r[1] = 0.0;
	r[2] = 0.0;
	
	wnd->Update(n, r, types, cell, rad);
	
	js = GetJSocket((char *)"localhost", port, 1);
	if( js == NULL ) {
		printf("Unavailable port, exiting OpenGL display system\n");
		exit(-1);
	}

	while (1) {
		connection = ListenJSocket(js, &caddr, &clen, 1, 0);
		if (connection > 0) {
			ReadJSocket(connection, cell, sizeof(double)*3 );
			ReadJSocket(connection, &new_n, sizeof(int) );
			// new n differs from old n?
			if (n != new_n) {
				free(types);
				free(r);
				free(rad);
				n = new_n;
				types = (int *)malloc(sizeof(int)*n);
				r = (double *)malloc(sizeof(double)*n*3);
				rad = (double *)malloc(sizeof(double)*n);
			}
			
			ReadJSocket(connection, rad, sizeof(double)*n);
			ReadJSocket(connection, types, sizeof(int)*n);
			ReadJSocket(connection, r, sizeof(double)*n*3);
			close(connection);
	
			wnd->Update(n, r, types, cell, rad);
		}
		if (wnd->ShouldQuitThread() == true) break;
	}
	
	FreeJSocket(js);
	
	free(types);
	free(r);
	free(rad);
	
	printf("Exiting program...\n");
	return NULL;
}


void *console_thread_func(void *data) {

	vis *v;
	char buffer[1024], *tokens[10];
	char *delims = (char *)" \t";
	int i, ntokens;
	char *p;
	
	v = (vis *)data;
	
	usleep(50000);
	printf("Molecular Simulation Visualization desktop\n");

	while (1) {
		printf("> ");
		p = fgets(buffer, 1024, stdin);
		if (p == NULL) continue;
		// strip eol
		buffer[strlen(buffer)-1] = '\0';
		
		ntokens = TokenizeString(buffer, delims, tokens, 10);
		if (ntokens < 1 || buffer[0] == '#' || tokens[0][0] == '#') continue;

		if ((strcasecmp(tokens[0],"h") == 0) || (strcasecmp(tokens[0], "?") == 0)) {
			printf("\nCommand menu:\n");
			printf("\t'h' or '?': shows this information screen\n");
			printf("\t'i': print all program information\n");
			printf("\t'v <x>': list visibility status of types, or toggle visibility of type x\n");
			printf("\t'q': end program. same as pressing 'q' in the visualisation window\n");
			printf("\n\nInside OpenGL window, use arrows up and down to zoom, left and right to rotate, 'q' to quit\n\n");
		}
		else if (strcasecmp(tokens[0],"v") == 0) {
			if (ntokens == 1) {
				printf("Non-visible types (otherwise visible):\n");
				v->BaseThreadLock();
					for (i = 0; i < MAX_TYPES; i++)
					{
						if (v->show_types[i] != 1) printf("%d\n", i);
					}
				v->BaseThreadUnlock();
			}
			else {
				i = atoi(tokens[1]);
				if (i < 0 || i >= MAX_TYPES) {
					printf("Bad type %d; allowed values are zero to %d\n", i, MAX_TYPES);
				}
				else {
					v->BaseThreadLock();
						v->show_types[i] = !v->show_types[i];
						v->data_invalidated = true;
						printf("Type %d has view status %d\n", i, v->show_types[i]);
					v->BaseThreadUnlock();
				}
			}
		}
		else if (strcasecmp(tokens[0],"i") == 0)
		{
			v->BaseThreadLock();
			for (i=0; i<MAX_TYPES; i++)
			{
				printf("Site type %d:\n", i);
				printf("\tambient: %f,%f,%f\n",  v->bead_amb[(i*4)+0],  v->bead_amb[(i*4)+1],  v->bead_amb[(i*4)+2]);
				printf("\tdiffuse: %f,%f,%f\n",  v->bead_diff[(i*4)+0], v->bead_diff[(i*4)+1], v->bead_diff[(i*4)+2]);
				printf("\tspecular: %f,%f,%f\n", v->bead_spec[(i*4)+0], v->bead_spec[(i*4)+1], v->bead_spec[(i*4)+2]);
				printf("\temit: %f,%f,%f\n",     v->bead_emit[(i*4)+0], v->bead_emit[(i*4)+1], v->bead_emit[(i*4)+2]);
				printf("\tshiny: %f\n", v->bead_shiny[i]);
				printf("\tshow: %d\n", v->show_types[i]);
			}
			v->BaseThreadUnlock();
		}
		else if (strcasecmp(tokens[0],"q") == 0) {
			exit( 0 );
		}
	}
	
	return NULL;
}

int main(int argc, char **argv)
{
//	int rc = 0;
	vis *v;
	jthread *console_thread;
	
	if (argc == 1) {
		printf("Usage: visualiser <listen_port_number>\n\n");
		exit(-1);
	}
	
	sscanf(argv[1], "%d", &port);

	v = new vis(20, argc, argv);

	console_thread = GetJThread(console_thread_func, v);
	if (console_thread == NULL) {
		printf("Unable to get create thread.\n");
		delete v;
		exit(-1);
	}
	RunJThread(console_thread);

// rc = v->Start(visFunc, 1000);
	v->Start(visFunc, 1000);
	FreeJThread(console_thread);
	
	return 0;
}
