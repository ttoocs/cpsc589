//Scott Saunders 10163541
// Based on the TA's template. (My own seems to be crashing as of writting)
#include <GLFW/glfw3.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <iostream>
#include <math.h>

#include <vector>
#include "Eigen/Dense"
#include "types.h"
using namespace std;

double zoom = 1;

GLFWwindow *window;
int w, h;
double mouseX, mouseY;
bool mouseDown;
bool mouseRDown;
bool shiftDown;

double visEps = 0.01;

std::vector<double> knots;
std::vector<vec2> points;
std::vector<double> weights;
int degree = 2;

int move_point = -1;

//#define INDFIX(X,Y) X[ std::max( (int) std::min(0,Y) , (int) X.size())]
#define INDFIX(X,Y) X[Y] 

int getLocPoint(double x, double y){
    int ret=-1;
    for(int i=0; i < points.size() ; i++){
      double len = (points[i] - vec2(mouseX,mouseY)).norm();
      if(len < visEps){
        ret=i;
        return ret;
      }
    }
  return -1;
}

void add_point(double x, double y, double weight=1){
  vec2 nvec = vec2(x,y);
  points.push_back(vec2(x,y));
  weights.push_back(weight);

//  std::cout << "x: " << x << "\ty: " << y << std::endl;
//  std::cout << "pnt:" <<std::endl;
//  std::cout << points[points.size()-1] << std::endl;
}

void del_point(){
  int i = getLocPoint(mouseX, mouseY);
  if( i != -1){
    points.erase(points.begin()+i);
    weights.erase(weights.begin()+i);
  }
}

void check_move(){
  //if(last_down = false && mouseRDown){
  if(move_point == -1){
    move_point = getLocPoint(mouseX,mouseY);
  }
  //last_down = mouseRDown;
}
void move_pnt(){
  if(move_point != -1 ){
    points[move_point] = vec2(mouseX, mouseY);
  }
}
void end_move(){
 // if(last_down && ! mouseDown && move_point != -1){
  if(move_point != -1 && mouseRDown){
    move_pnt();
    move_point = -1;
  }
}

double def_color[3] = {1,1,1};
void draw_circle(double r, double xoff, double yoff, double color[3]){
  bool first = true;
  double lastPosX;
  double lastPosY;
  double stepSize = 1.0/10.0;
  for(double u=0; u< 2*3.1415926535; u+=stepSize){
    double X=sin(u)*r;
    double Y=cos(u)*r;
    X +=xoff;
    Y +=yoff;

    if(! first){
      glBegin(GL_LINES);
      glColor3f(color[0],color[1],color[2]);
      glVertex2f(lastPosX,lastPosY);
      glVertex2f(X,Y);
      glEnd();
    }
    lastPosX=X;
    lastPosY=Y;
    first=false;
  }
}


std::vector<double> gen_stdr(){
  int len = points.size();
  std::vector<double> knots;

  if(points.size() < 2)
    return(knots);

  for(int i = 0; i < degree  ; i++){
    knots.push_back(0);
  }

  double stepSize = 1 / ( (double) ( (double) (points.size() -1) - (double) (degree + 1) +2)) ;
//  std::cout << points.size() << "," << degree << "," << std::endl;
//  std::cout << "stepSize = " << stepSize <<"\t inv: " << 1/stepSize << std::endl;
  for(double i = 0; i < 1 ; i+=stepSize ){
    knots.push_back(i);
  }

  for(int i = 0; i < degree  ; i++){
    knots.push_back(1);
  }
  return knots;
}


int getDelta(double u){
  if(knots.size() < 2)
    return -1;

  for(int i=0; i < knots.size() -1 ; i++){
    //if(u >= knots[i] && ( ( u < knots[i+1] ) || (  ((i+1) >= knots.size()) && ((i+1) > knots[knots.size()]) )))
    if(u >= INDFIX(knots,i) && u < INDFIX(knots,i+1)) // || (  ((i+1) >= knots.size()) && ((i+1) > knots[knots.size()]) )))
      return i;
  }
  if(u > knots[knots.size()-1])
    return knots.size();
  std::cout << "Delta broke" << std::endl;
  return -1;
}

//There is some error occuring in here, about an edge case,
//It's noticable when the degree =1.
double bruteN(int i, int r, double u){
  if( r < 1)
    std::cout << "zero r" << std::endl;

  if ( r == 1){
  int s = knots.size();
  if ( i >= 0 && s > i  && INDFIX(knots,i) <= u ){

    if( ( u < INDFIX(knots,i+1) )){ // && ((i+1) < knots.size()) && (i+1 > 0 ) ){
      return 1;
    }
  }

//    std::cout << " if ( u_i \\leq  u \\leq u_{i+1} ): 1 \t";// << std::endl;
    return 0;
  } 
  
 
//  std::cout << " ( "  ;
//  std::cout << " ( \\frac { u - u_{" << i << "} }{ u_{ " << i + r - 1 << "} - u_{" << i << "} }) " ;
  double LHN = (( u- INDFIX(knots,i) ));
  double LHD = ( INDFIX(knots,i+r-1) - INDFIX(knots,i));
  double LR = bruteN(i,r-1,u);
  double  L = (LHN/LHD) * LR;

//  std::cout << "\n + " ;
//  std::cout << " ( \\frac { u_{" << i+r << "} - u }{ u_{ " << i + r  << "} - u_{" << i+1 << "} }) " ;
  double RHN = (( INDFIX(knots,i+r) - u)) ;
  double RHD =(INDFIX(knots,i+r) - INDFIX(knots,i+1));
  double RR = bruteN(i+1,r-1,u);
  double R= (RHN/RHD)*RR;

//  std::cout << " ) ";
  if( LHD == 0  ||RHD == 0 )
    return 0;
  if( (i+1) > knots.size())
    return 0;
  
//  std::cout << "LHD:\t" << L << "\t\tRHD:\t" << R << std::endl;
  return L+R;
}

vec2 bruteS(double u){
  vec2 s= vec2(0,0);
  //Here we can add in an optimzatoin
  //for(int i = min(((int) u) - degree,0) 
  for(int i = 0; i < points.size(); i++){
    s=s + points[i]*bruteN(i,degree,u);
  }
  return s;
}
vec2 fancy(double u){
    int delta = getDelta(u);
//    int delta = (int) u;
    std::vector<vec2> coeff;

    for(int i=0 ; i < degree; i++){
      coeff.push_back(points[delta-1]);
    }

    for(int r=degree+1; r == 2 ; r -= 1){
      int i = delta;
      for(int s=0 ; s == r-2; s++){
        double omega = (u - INDFIX(knots,i))/(INDFIX(knots,i+r-1) - INDFIX(knots,i));
        coeff[s]  = omega * coeff[s] + (1-omega)*coeff[s+1];
        i=i-1;
      }
    }
    return coeff[0];
}

std::vector<vec2> get_NURB(){ 
  //NOTE: Not yet effishent
  //NOTE: Not yet NURBS

  double stepSize = 0.001;

  knots = gen_stdr();
  std::vector<vec2> pos;

  if(points.size() < degree)
    return pos;

  
  for(double u=0 ; u < 1  ; u+= stepSize){
    pos.push_back(bruteS(u));
//    pos.push_back(fancy(u));
  }
  
  

  return pos;
}


void render() {
	glEnable (GL_DEPTH_TEST);
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Functions for changing transformation matrix
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity ();
	glTranslatef(0.0f, 0.0f, 0.0f );
	glRotatef (0, 0.0f, 0.0f, 1.0f);
	glScalef (zoom,zoom, zoom);

	//Functions for changing projection matrix
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho (-1, 1, -1, 1, -1, 1);
	//gluPerspective (fov, aspect ratio, near plane, far plane)
	//glFrustum


  move_pnt();
  for(int i = 0; i < points.size() ; i++){
    draw_circle(visEps*weights[i], points[i](0), points[i](1),def_color);
  }

  std::vector<vec2> curve;
  curve = get_NURB();
  //glBegin(GL_LINE_STRIP); 
  for(int i = 1; i < curve.size(); i++){
    glBegin(GL_LINES);
    glVertex2f(curve[i-1](0), curve[i-1](1));
    glVertex2f(curve[i](0), curve[i](1));
    glEnd();
  }
  //glEnd();

}

void keyboard (GLFWwindow *sender, int key, int scancode, int action, int mods) {

  if (action == GLFW_RELEASE){
//    if(mods == GLFW_MOD_SHIFT)
    if(key == GLFW_KEY_LEFT_SHIFT || key == GLFW_KEY_RIGHT_SHIFT)
      shiftDown = false;

    return;
  }else if (key == GLFW_KEY_Q){
		cout << "Q: Quitting! " << std::endl;
    glfwSetWindowShouldClose(window, GLFW_TRUE);
  }

  if(key == GLFW_KEY_A && action == GLFW_PRESS){
    degree++;
    std::cout << "Degree set to: " << degree << std::endl;
  }
  if(key == GLFW_KEY_Z && action == GLFW_PRESS){
    degree--;
    if(degree < 1)
      degree = 1;
    std::cout << "Degree set to: " << degree << std::endl;
  }
  

//  if(mods & GLFW_MOD_SHIFT != 0){
  if(key == GLFW_KEY_LEFT_SHIFT || key == GLFW_KEY_RIGHT_SHIFT){
    shiftDown = true;
//    std::cout << "SHIFT" << std::endl;
  }
}

void mouseClick (GLFWwindow *sender, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS){
//		cout << mouseX << ' ' << mouseY << '\n';
    mouseDown = true;
    add_point(mouseX, mouseY);
  }
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE){
    mouseDown = false;
    end_move();
  }
	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS){
//		cout << mouseX << ' ' << mouseY << '\n';
    if(shiftDown){
      del_point();
    }else{
      mouseRDown = true;
      check_move();
    }
  }
  if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE){
    mouseDown = false;
    end_move();
  }


}
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset){
  int i = getLocPoint(mouseX,mouseY);
  if( i != -1)
    weights[i] += xoffset/50;

}


void mousePos (GLFWwindow *sender, double x, double y) {
	mouseX = (x/(w/2.0) - 1);
	mouseY = 1 - (y/(h/2.0));
  //Localized to mouse-pose.
}


int main () {
	if (!glfwInit())
		return 1;

	window = glfwCreateWindow (640, 640, "Scott Saunders A2", NULL, NULL);
	if (!window)
		return 1;

	glfwMakeContextCurrent (window);
	glfwSetKeyCallback (window, keyboard);
	glfwSetMouseButtonCallback (window, mouseClick);
	glfwSetCursorPosCallback (window, mousePos);
  glfwSetScrollCallback(window, scroll_callback);
 
  points.clear();
  weights.clear();
  knots.clear();
  knots = gen_stdr();
//  bruteN(0,3,2);  

//  for(int i = 0; i < knots.size(); i++){
//    std::cout << "knot[ " << i << "]: " << knots[i] << std::endl;
//  }
	while (!glfwWindowShouldClose (window)) {
		glfwGetFramebufferSize (window, &w, &h);
		glViewport (0, 0, w, h);

		render();

		glfwSwapBuffers (window);
//    if(animate)
//      glfwPollEvents();
//    else
  		glfwWaitEvents(); //No need to animate if not animating.
	}

	glfwDestroyWindow (window);
	glfwTerminate();
	return 0;
}


