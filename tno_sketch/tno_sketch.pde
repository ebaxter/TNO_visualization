String data_dir = "../data/";

int npart; //number of TNOs
int nplanet = 9; //number of planets

float[] xx; //positions of objects
float[] yy; 
float[] zz;
float[] xdot; //velocities of objects
float[] ydot; 
float[] zdot;
float[] r_arr; //radii (diameter?) of objects [km]
float[] a_arr; //semi major axis [m]
float[] e_arr; //eccentricity
float[] lilomega_arr; //argument of periapsis [radians]
float[] bigomega_arr; //longitude of ascending node [radians]
float[] inc_arr; //inclination [radians]
float[] M0_arr; //Mean anomaly
float t;
float dt;
float t0;

//Scaling of plot relative to physical units
float scaling;

//Camera position
float centerx;  
float centery;
float centerz;

//For keeping track of mouse position
boolean tracking;
float mousestartx = 0.;
float mousestarty = 0.;
float mx;
float my;

//Useful constants
float pi = 3.14159265;  //Constant
float mu = 1.32712440041e20; 
float AU_to_m = 1.49597870691e11;  // Unit conversion
float sun_radius = 5.0;  //how large the Sun appears

void setup() {
  size(700, 700, P3D);
  
  //Simulation settings
  scaling = 6.0;  //How to scale orbital distances
  float radius_scaling = 150.0;  //How large to make the spheres when plotting
  dt = 10.0*0.3*10000.0/radius_scaling;
  float minrad = 30.0; //Minimum size of TNO to include [km]
  
  centerx = 0.0;
  centery = 0.0;
  centerz = 0.0;
  tracking = false;
  
  //initial conditions
  t0 = 0.0; //Julian day
  t = 0.0; //Julian day

  //Data files
  String[] lines = loadStrings(data_dir + "Distant.txt");
  String[] rad_lines = loadStrings(data_dir + "tno_radii.dat");
  String[] planet_lines = loadStrings(data_dir + "planet_orbital.dat");

  //number of TNOs
  npart = lines.length - 1;

  //Get orbital elements and radii
  int nstart = 500;
  float a = 0.0;
  float e = 0.0;
  float a_au = 0.0;
  float rtemp = 0.0;
  float lilomega = 0.0;
  float bigomega = 0.0;
  float M0 = 0.0;
  float inc = 0.0;
  float id;
 
  float[] a_arr_temp = new float[nstart]; 
  float[] e_arr_temp = new float[nstart]; 
  float[] lilomega_arr_temp = new float[nstart]; 
  float[] bigomega_arr_temp = new float[nstart]; 
  float[] inc_arr_temp = new float[nstart]; 
  float[] M0_arr_temp = new float[nstart]; 
  float[] r_arr_temp = new float[nstart];

  //Load TNO data
  float maxrad = 0.0;
  float maxa = 0.0;
  int i = 0;
  for (int ii = 0; ii < lines.length; ii++) {
    id = float(lines[ii].substring(1-1,7-1));
    rtemp = 20.0; //default value
    for (int j = 0; j < rad_lines.length; j++){
        String[] rad_split = split(rad_lines[j]," ");
        if (float(rad_split[0]) ==  id){
            rtemp = float(rad_split[1]); 
        }
    }
    a_au = float(lines[ii].substring(93-1,103-1));
    e = float(lines[ii].substring(71-1,79-1));
    lilomega = float(lines[ii].substring(38-1,46-1))*pi/180.0;
    bigomega = float(lines[ii].substring(49-1,57-1))*pi/180.0;
    inc = float(lines[ii].substring(60-1,68-1))*pi/180.0;
    M0 = float(lines[ii].substring(27-1,35-1))*pi/180.0;
    
    //For testing purposes
    boolean test = false;
    if (test){
       a_au = 10.0 + 0.1*(i+1.0);
       e = 0.0;
       lilomega = 0.0;
       bigomega = 0.0;
       inc = 0.0;
       M0 = 0.0;
    }
    
    //Only include the TNO if radius is above minimum
    if (rtemp > minrad & i < npart){
      a = a_au*AU_to_m;
      
      r_arr_temp[i] = rtemp/radius_scaling;
      a_arr_temp[i] = a;
      e_arr_temp[i] = e;
      lilomega_arr_temp[i] = lilomega;
      bigomega_arr_temp[i] = bigomega;
      inc_arr_temp[i] = inc;
      M0_arr_temp[i] = M0;
      i++;
    }

    //Keep track of minimum and maximum radii
    if (rtemp > maxrad){
     maxrad = rtemp; 
    }
    if (a_au > maxa){
      maxa = a_au;
    }
  }
  println("maxrad = ", maxrad);
  println("maxa = ", maxa);
  println("Number of bodies above min rad = ", i);
  npart = i;
 
  //Load planet data
  //http://www.met.rdg.ac.uk/~ross/Astronomy/Planets.html
  float[] a_arr_pl = new float[nplanet]; //semi major axis, m
  float[] e_arr_pl = new float[nplanet]; //eccentricity
  float[] lilomega_arr_pl = new float[nplanet]; //argument of periapsis [radians]
  float[] bigomega_arr_pl = new float[nplanet]; //longitude of ascending node [radians]
  float[] inc_arr_pl = new float[nplanet]; //inclination [radians]
  float[] M0_arr_pl = new float[nplanet]; //Mean anomaly
  float[] r_arr_pl = new float[nplanet]; //Diameter [m?]
  
  for (int jj = 1; jj <= nplanet; jj++){
      String[] pl_split = split(planet_lines[jj],",");
      a_arr_pl[jj-1] = float(pl_split[1])*AU_to_m;
      e_arr_pl[jj-1] = float(pl_split[2]);
      inc_arr_pl[jj-1] = float(pl_split[3])*pi/180.;
      bigomega_arr_pl[jj-1] = float(pl_split[3])*pi/180.;
      lilomega_arr_pl[jj-1] = float(pl_split[4])*pi/180.;
      //Is this correct?  Is MO same as mean longitude?
      M0_arr_pl[jj-1] = float(pl_split[5])*pi/180.;
      r_arr_pl[jj-1] = 2.0;
  }

  //Combine TNO and planet data into single arrays
  a_arr = new float[npart + nplanet]; //semi major axis, m
  e_arr = new float[npart + nplanet]; //eccentricity
  lilomega_arr = new float[npart + nplanet]; //argument of periapsis [radians]
  bigomega_arr = new float[npart + nplanet]; //longitude of ascending node [radians]
  inc_arr = new float[npart + nplanet]; //inclination [radians]
  M0_arr = new float[npart + nplanet]; //Mean anomaly
  r_arr = new float[npart + nplanet]; //Diameter [m?]
  for (int j = 0; j < npart + nplanet; j++){
      if (j < nplanet){
        r_arr[j] = r_arr_temp[j];
        a_arr[j] = a_arr_pl[j];
        e_arr[j] = e_arr_pl[j];
        lilomega_arr[j] = lilomega_arr_pl[j];
        bigomega_arr[j] = bigomega_arr_pl[j];
        inc_arr[j] = inc_arr_pl[j];
        M0_arr[j] = M0_arr_pl[j];
        r_arr[j] = r_arr_pl[j];
      }
      if (j >= nplanet){
        r_arr[j] = r_arr_temp[j-nplanet];
        a_arr[j] = a_arr_temp[j-nplanet];
        e_arr[j] = e_arr_temp[j-nplanet];
        lilomega_arr[j] = lilomega_arr_temp[j-nplanet];
        bigomega_arr[j] = bigomega_arr_temp[j-nplanet];
        inc_arr[j] = inc_arr_temp[j-nplanet];
        M0_arr[j] = M0_arr_temp[j-nplanet];
      }
  }

  //initialize positions
  xx = new float[npart + nplanet];
  yy = new float[npart + nplanet];
  zz = new float[npart + nplanet];
  xdot = new float[npart + nplanet];
  ydot = new float[npart + nplanet];
  zdot = new float[npart + nplanet];
  update_positions(t, t0);

  frameRate(20);
  noStroke();  
  camera(width/2.0, height/2.0, -500.0, 0., 0., 0., 0, 1, 0);
}

void update_positions(float t, float t0) {
  float[] rt = new float[3];
  float[] rtdot = new float[3];

  for (int i = 0; i < npart + nplanet; i++) {
    rt = get_pos(t, t0, M0_arr[i], a_arr[i], e_arr[i], lilomega_arr[i], bigomega_arr[i], inc_arr[i]);
    //rtdot = get_posdot(t, t0, M0_arr[i], a_arr[i], e_arr[i], lilomega_arr[i], bigomega_arr[i], inc_arr[i]);
    xx[i] = rt[0];
    yy[i] = rt[1];
    zz[i] = rt[2];
    //xdot[i] = rtdot[0];
    //ydot[i] = rtdot[1];
    //zdot[i] = rtdot[2];
  }
}

float[] get_pos(float t, float t0, float M0, float a, float e, float lilomega, float bigomega, float inclination) {
  float Mt = M_of_t(t, a, t0, M0);
  float Et = E_of_t(Mt, e);
  float nut = nu_of_t(e, Et);
  float rct = rc_of_t(a, e, Et); 
  float[] ot = get_ot(rct, nut);
  float[] rt = get_rt_au(ot, bigomega, inclination, lilomega, t);

  return rt;
}

float[] get_posdot(float t, float t0, float M0, float a, float e, float lilomega, float bigomega, float inclination) {  
  float Mt = M_of_t(t, a, t0, M0);
  float Et = E_of_t(Mt, e);
  float nut = nu_of_t(e, Et);
  float rct = rc_of_t(a, e, Et);
  float[] odot = get_otdot(a, rct, Et, e);
  float[] rtdot = get_rtdot_au(odot, bigomega, inclination, lilomega, t);

  return rtdot;
}

float M_of_t(float t, float a, float t0, float M0) {
  float M;
  if (t == t0) {
    M = M0;
  } else {
    float Deltat = 86400.*(t - t0); 
    //In this form to prevent rounding errors
    M = M0 + Deltat*pow(AU_to_m, -3./2.)*sqrt(mu/pow(a/AU_to_m, 3));
  }

  //Normalize
  if (M > 2.*pi) {
    M = M - 2.*pi*floor(M/(2.*pi));
  }
  return M;
}

float E_of_t(float Mt, float e)
{
  float Etprevious = 0.0;
  float Et = Mt;
  float desired_precision = 0.00001;
  float current_precision = 1.0;
  float newterm = 0.0;
  int num_iter = 0;
  while (current_precision > desired_precision) {
    Etprevious = Et;
    newterm = (Et - e*sin(Et) - Mt)/(1.0 - e*cos(Et));
    Et = Et - newterm;
    num_iter++;
    current_precision = (Et - Etprevious)/Etprevious;
  }
  //println("num iter = ", num_iter);
  return Et;
}

float nu_of_t(float e, float Et)
{
  float nu = 2.*arctan2(sqrt(1.0 + e)*sin(Et/2.0), sqrt(1-e)*cos(Et/2.0));
  return nu;
}

float rc_of_t(float a, float e, float Et) {
  return a*(1.0-e*cos(Et));
}

float[] get_ot(float rc, float nut) {
  float[] ot = new float[3];
  ot[0] = rc*cos(nut);
  ot[1] = rc*sin(nut);
  ot[2] = 0.0;
  return ot;
}

float[] get_otdot(float a, float rc, float Et, float e) {
  float[] otdot = new float[3];
  float prefactor = sqrt(mu*a)/rc;
  otdot[0] = prefactor*(-sin(Et));
  otdot[1] = prefactor*(sqrt(1.0-e*e)*cos(Et));
  otdot[2] = 0.0;
  return otdot;
}

float[] get_rt_au(float[] ovec, float bigO, float inc, float lilo, float t) {
  float[] rt = new float[3];
  float[] rt_au = new float[3];

  rt[0] = ovec[0]*(cos(lilo)*cos(bigO) - sin(lilo)*cos(inc)*sin(bigO)) - ovec[1]*(sin(lilo)*cos(bigO) + cos(lilo)*cos(inc)*sin(bigO));
  rt[1] = ovec[0]*(cos(lilo)*sin(bigO) + sin(lilo)*cos(inc)*cos(bigO)) + ovec[1]*(cos(lilo)*cos(inc)*cos(bigO) - sin(lilo)*sin(bigO));
  rt[2] = ovec[0]*(sin(lilo)*sin(inc)) + ovec[1]*(cos(lilo)*sin(inc));

  rt_au[0] = rt[0]/AU_to_m;
  rt_au[1] = rt[1]/AU_to_m;
  rt_au[2] = rt[2]/AU_to_m;
  return rt_au;
}

float[] get_rtdot_au(float[] odotvec, float bigO, float inc, float lilo, float t) {
  float[] rtdot = new float[3];
  float[] rtdot_au = new float[3];

  rtdot[0] = odotvec[0]*(cos(lilo)*cos(bigO) - sin(lilo)*cos(inc)*sin(bigO)) - odotvec[1]*(sin(lilo)*cos(bigO) + cos(lilo)*cos(inc)*sin(bigO));
  rtdot[1] = odotvec[0]*(cos(lilo)*sin(bigO) + sin(lilo)*cos(inc)*cos(bigO)) + odotvec[1]*(cos(lilo)*cos(inc)*cos(bigO) - sin(lilo)*sin(bigO));
  rtdot[2] = odotvec[0]*(sin(lilo)*sin(inc)) + odotvec[1]*(cos(lilo)*sin(inc));

  rtdot_au[0] = rtdot[0]/AU_to_m;
  rtdot_au[1] = rtdot[1]/AU_to_m;
  rtdot_au[2] = rtdot[2]/AU_to_m;
  return rtdot_au;
}


void draw() {
  background(0);
  lights();
  directionalLight(128, 128, 128, 0, 0, 1);

  //Determine camera position and angle, etc.
  if (!tracking) {
    mx = mouseX;
    my = mouseY;
  }

  //fixed view
  boolean do_moving_view = true;
  if (!do_moving_view) {
    camera(width/2.0, height/2.0, -500.0, 0., 0., 0., 0, 1, 0);
    //rotateY( - PI/4 + 0.2);
  }

  //Moving view
  if (do_moving_view) {
    float cameraY = height/2.0;
    float fov = my/float(width) * PI/2;
    float cameraZ = cameraY / tan(fov / 2.0);
    float aspect = float(width)/float(height);
    sun_radius = fov*10.0;

    rotateY(mx/float(width)* PI - PI/4 + 0.2);  //This setting gets nice agreement between mouse position and viewing angle
    perspective(fov, aspect, cameraZ/10.0, cameraZ*10.0);
  }

  if (mousePressed) {
    if (!tracking) {
      mousestartx = mouseX;
      mousestarty = mouseY;
      tracking = true;
    }
    centerx = mouseX - mousestartx;
    centery = mouseY - mousestarty;
  }

  //Draw the sun
  translate(0., 0., 0.);
  translate(centerx, centery);
  color sun_color = #FFE600;
  fill(sun_color);
  sphere(sun_radius); 
  color planet_color = #F44242;
  fill(planet_color);
 
  //Draw the spheres
  for (int i = 0; i < npart + nplanet; i++) {
    if (i >= nplanet){
      color tno_color = #00BCFF;
      fill(tno_color);
    }
    translate(xx[i]*scaling, yy[i]*scaling, zz[i]*scaling);
    //Draw object
    sphere(r_arr[i]);   
    //Undo translation
    translate(-xx[i]*scaling, -yy[i]*scaling, -zz[i]*scaling);
  }

  t = t + dt;
  update_positions(t, t0);
}

void mouseReleased() {
  tracking = false;
}

float arctan2(float y, float x) {
  float returnval = 0.0  ;
  if (x > 0) {
    returnval = atan(y/x);
  }
  if (y >= 0 && x < 0) {
    returnval = atan(y/x) + pi;
  }
  if (y < 0 && x < 0) {
    returnval = atan(y/x) - pi;
  }
  if (y > 0 && x == 0) {
    returnval = pi/2.0;
  }
  if (y < 0 && x == 0) {
    returnval = -pi/2.0;
  }
  return returnval;
}
