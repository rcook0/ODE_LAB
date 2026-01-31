
use wasm_bindgen::prelude::*;
use ode_solvers::*;
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize)]
pub struct Series {
    pub t: Vec<f64>,
    pub y: Vec<f64>, // flattened
    pub dim: usize,
}


type RhsFn = fn(t: f64, y: &[f64], dy: &mut [f64], params: &std::collections::BTreeMap<String, f64>);

fn integrate_rk4_vec(
    rhs: RhsFn,
    params: &std::collections::BTreeMap<String, f64>,
    mut t: f64,
    t1: f64,
    mut y: Vec<f64>,
    dt: f64,
) -> Series {
    let n = y.len();
    let mut ts: Vec<f64> = Vec::new();
    let mut ys: Vec<f64> = Vec::new();

    let mut k1 = vec![0.0; n];
    let mut k2 = vec![0.0; n];
    let mut k3 = vec![0.0; n];
    let mut k4 = vec![0.0; n];
    let mut yt = vec![0.0; n];

    let mut h = if dt == 0.0 { 1e-3 } else { dt };
    let max_steps = ((t1 - t).abs() / h.abs()).ceil() as usize + 2;

    ts.push(t);
    ys.extend_from_slice(&y);

    for _ in 0..max_steps {
        // stop if we've reached or passed t1 in the direction of travel
        if (h > 0.0 && t >= t1) || (h < 0.0 && t <= t1) { break; }
        if (t + h - t1).signum() == h.signum() && (t + h - t1).abs() > 1e-12 {
            h = t1 - t; // final step
        }

        rhs(t, &y, &mut k1, params);

        for i in 0..n { yt[i] = y[i] + 0.5*h*k1[i]; }
        rhs(t + 0.5*h, &yt, &mut k2, params);

        for i in 0..n { yt[i] = y[i] + 0.5*h*k2[i]; }
        rhs(t + 0.5*h, &yt, &mut k3, params);

        for i in 0..n { yt[i] = y[i] + h*k3[i]; }
        rhs(t + h, &yt, &mut k4, params);

        for i in 0..n {
            y[i] = y[i] + (h/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        }
        t += h;

        ts.push(t);
        ys.extend_from_slice(&y);
    }

    Series { t: ts, y: ys, dim: n }
}

// -------- Per-model RHS functions (Vec-based, dimension N) --------

fn rhs_linear_2x2(_t: f64, y: &[f64], dy: &mut [f64], p: &std::collections::BTreeMap<String, f64>) {
    let a11 = *p.get("a11").unwrap_or(&0.0);
    let a12 = *p.get("a12").unwrap_or(&0.0);
    let a21 = *p.get("a21").unwrap_or(&0.0);
    let a22 = *p.get("a22").unwrap_or(&0.0);
    let x = y[0]; let yy = y[1];
    dy[0] = a11*x + a12*yy;
    dy[1] = a21*x + a22*yy;
}

fn rhs_vanderpol(_t: f64, y: &[f64], dy: &mut [f64], p: &std::collections::BTreeMap<String, f64>) {
    let mu = *p.get("mu").unwrap_or(&1.0);
    let x = y[0]; let v = y[1];
    dy[0] = v;
    dy[1] = mu*(1.0 - x*x)*v - x;
}

fn rhs_predator_prey(_t: f64, y: &[f64], dy: &mut [f64], p: &std::collections::BTreeMap<String, f64>) {
    let alpha = *p.get("alpha").unwrap_or(&1.0);
    let beta  = *p.get("beta").unwrap_or(&0.1);
    let delta = *p.get("delta").unwrap_or(&0.075);
    let gamma = *p.get("gamma").unwrap_or(&1.5);
    let x = y[0]; let yy = y[1];
    dy[0] = alpha*x - beta*x*yy;
    dy[1] = delta*x*yy - gamma*yy;
}

fn rhs_duffing(t: f64, y: &[f64], dy: &mut [f64], p: &std::collections::BTreeMap<String, f64>) {
    let delta = *p.get("delta").unwrap_or(&0.2);
    let alpha = *p.get("alpha").unwrap_or(&-1.0);
    let beta  = *p.get("beta").unwrap_or(&1.0);
    let gamma = *p.get("gamma").unwrap_or(&0.3);
    let omega = *p.get("omega").unwrap_or(&1.2);
    let x = y[0]; let v = y[1];
    dy[0] = v;
    dy[1] = -delta*v - alpha*x - beta*x*x*x + gamma*(omega*t).cos();
}

fn rhs_sir(_t: f64, y: &[f64], dy: &mut [f64], p: &std::collections::BTreeMap<String, f64>) {
    let beta = *p.get("beta").unwrap_or(&0.5);
    let gamma = *p.get("gamma").unwrap_or(&0.2);
    let s = y[0]; let i = y[1];
    dy[0] = -beta*s*i;
    dy[1] = beta*s*i - gamma*i;
    dy[2] = gamma*i;
}

fn rhs_lorenz63(_t: f64, y: &[f64], dy: &mut [f64], p: &std::collections::BTreeMap<String, f64>) {
    let sigma = *p.get("sigma").unwrap_or(&10.0);
    let rho = *p.get("rho").unwrap_or(&28.0);
    let beta = *p.get("beta").unwrap_or(&(8.0/3.0));
    let x=y[0]; let yy=y[1]; let z=y[2];
    dy[0] = sigma*(yy-x);
    dy[1] = x*(rho-z) - yy;
    dy[2] = x*yy - beta*z;
}

fn rhs_rossler(_t: f64, y: &[f64], dy: &mut [f64], p: &std::collections::BTreeMap<String, f64>) {
    let a = *p.get("a").unwrap_or(&0.2);
    let b = *p.get("b").unwrap_or(&0.2);
    let c = *p.get("c").unwrap_or(&5.7);
    let x=y[0]; let yy=y[1]; let z=y[2];
    dy[0] = -yy - z;
    dy[1] = x + a*yy;
    dy[2] = b + z*(x - c);
}

fn rhs_chua(_t: f64, y: &[f64], dy: &mut [f64], p: &std::collections::BTreeMap<String, f64>) {
    let alpha = *p.get("alpha").unwrap_or(&15.6);
    let beta  = *p.get("beta").unwrap_or(&28.0);
    let m0    = *p.get("m0").unwrap_or(&-1.143);
    let m1    = *p.get("m1").unwrap_or(&-0.714);
    let x=y[0]; let yy=y[1]; let z=y[2];
    let fx = m1*x + 0.5*(m0-m1)*((x+1.0).abs() - (x-1.0).abs());
    dy[0] = alpha*(yy - x - fx);
    dy[1] = x - yy + z;
    dy[2] = -beta*yy;
}

fn rhs_thomas(_t: f64, y: &[f64], dy: &mut [f64], p: &std::collections::BTreeMap<String, f64>) {
    let a = *p.get("a").unwrap_or(&0.208186);
    let x=y[0]; let yy=y[1]; let z=y[2];
    dy[0] = yy.sin() - a*x;
    dy[1] = z.sin()  - a*yy;
    dy[2] = x.sin()  - a*z;
}

fn params_from_js(params: JsValue) -> std::collections::BTreeMap<String, f64> {
    if params.is_undefined() || params.is_null() {
        return std::collections::BTreeMap::new();
    }
    serde_wasm_bindgen::from_value(params).unwrap_or_else(|_| std::collections::BTreeMap::new())
}

fn integrate_model_vec(model: &str, params: &std::collections::BTreeMap<String, f64>, y0: Vec<f64>, t0: f64, t1: f64, dt: f64) -> Series {
    match model {
        "linear" => integrate_rk4_vec(rhs_linear_2x2, params, t0, t1, y0, dt),
        "vdp" => integrate_rk4_vec(rhs_vanderpol, params, t0, t1, y0, dt),
        "lv" => integrate_rk4_vec(rhs_predator_prey, params, t0, t1, y0, dt),
        "duffing" => integrate_rk4_vec(rhs_duffing, params, t0, t1, y0, dt),
        "sir" => integrate_rk4_vec(rhs_sir, params, t0, t1, y0, dt),
        "lorenz" => integrate_rk4_vec(rhs_lorenz63, params, t0, t1, y0, dt),
        "rossler" => integrate_rk4_vec(rhs_rossler, params, t0, t1, y0, dt),
        "chua" => integrate_rk4_vec(rhs_chua, params, t0, t1, y0, dt),
        "thomas" => integrate_rk4_vec(rhs_thomas, params, t0, t1, y0, dt),
        _ => integrate_rk4_vec(rhs_linear_2x2, params, t0, t1, y0, dt),
    }
}

#[wasm_bindgen]
pub fn integrate_model(model: String, params: JsValue, y0: Vec<f64>, t0: f64, t1: f64, dt: f64) -> JsValue {
    let p = params_from_js(params);
    let series = integrate_model_vec(&model, &p, y0, t0, t1, dt);
    serde_wasm_bindgen::to_value(&series).unwrap()
}

fn flatten<const N: usize>(ys: &Vec<[f64; N]>) -> Vec<f64> {
    ys.iter().flat_map(|v| v.iter().cloned()).collect()
}

struct Lin2 { a11:f64,a12:f64,a21:f64,a22:f64 }
impl System<[f64;2]> for Lin2 {
    fn system(&self,_t:f64,y:&[f64;2],dy:&mut[f64;2]){
        dy[0]=self.a11*y[0]+self.a12*y[1];
        dy[1]=self.a21*y[0]+self.a22*y[1];
    }
}
#[wasm_bindgen]
pub fn integrate_linear_2x2(a11:f64,a12:f64,a21:f64,a22:f64,x0:f64,y0:f64,t0:f64,t1:f64,dt:f64)->JsValue{
    let mut p = std::collections::BTreeMap::new();
    p.insert("a11".into(), a11); p.insert("a12".into(), a12);
    p.insert("a21".into(), a21); p.insert("a22".into(), a22);
    let series = integrate_model_vec("linear", &p, vec![x0,y0], t0, t1, dt);
    serde_wasm_bindgen::to_value(&series).unwrap()
}


struct Vdp{ mu:f64 }
impl System<[f64;2]> for Vdp {
    fn system(&self,_t:f64,y:&[f64;2],dy:&mut[f64;2]){
        let x=y[0]; let v=y[1];
        dy[0]=v;
        dy[1]=self.mu*(1.0-x*x)*v - x;
    }
}
#[wasm_bindgen]
pub fn integrate_vanderpol(mu:f64,x0:f64,y0:f64,t0:f64,t1:f64,dt:f64)->JsValue{
    let mut p = std::collections::BTreeMap::new();
    p.insert("mu".into(), mu);
    let series = integrate_model_vec("vdp", &p, vec![x0,y0], t0, t1, dt);
    serde_wasm_bindgen::to_value(&series).unwrap()
}


struct LV{ alpha:f64,beta:f64,delta:f64,gamma:f64 }
impl System<[f64;2]> for LV {
    fn system(&self,_t:f64,y:&[f64;2],dy:&mut[f64;2]){
        let x=y[0]; let p=y[1];
        dy[0]=self.alpha*x - self.beta*x*p;
        dy[1]=self.delta*x*p - self.gamma*p;
    }
}
#[wasm_bindgen]
pub fn integrate_predator_prey(alpha:f64,beta:f64,delta:f64,gamma:f64,x0:f64,y0:f64,t0:f64,t1:f64,dt:f64)->JsValue{
    let mut p = std::collections::BTreeMap::new();
    p.insert("alpha".into(), alpha); p.insert("beta".into(), beta);
    p.insert("delta".into(), delta); p.insert("gamma".into(), gamma);
    let series = integrate_model_vec("lv", &p, vec![x0,y0], t0, t1, dt);
    serde_wasm_bindgen::to_value(&series).unwrap()
}


struct Duff{ delta:f64,alpha:f64,beta:f64,gamma:f64,omega:f64 }
impl System<[f64;2]> for Duff {
    fn system(&self,t:f64,y:&[f64;2],dy:&mut[f64;2]){
        let x=y[0]; let v=y[1];
        dy[0]=v;
        dy[1]=-self.delta*v - self.alpha*x - self.beta*x*x*x + self.gamma*(self.omega*t).cos();
    }
}
#[wasm_bindgen]
pub fn integrate_duffing(delta:f64,alpha:f64,beta:f64,gamma:f64,omega:f64,x0:f64,y0:f64,t0:f64,t1:f64,dt:f64)->JsValue{
    let mut p = std::collections::BTreeMap::new();
    p.insert("delta".into(), delta); p.insert("alpha".into(), alpha);
    p.insert("beta".into(), beta); p.insert("gamma".into(), gamma); p.insert("omega".into(), omega);
    let series = integrate_model_vec("duffing", &p, vec![x0,y0], t0, t1, dt);
    serde_wasm_bindgen::to_value(&series).unwrap()
}


struct Sir{ beta:f64,gamma:f64 }
impl System<[f64;3]> for Sir {
    fn system(&self,_t:f64,y:&[f64;3],dy:&mut[f64;3]){
        let s=y[0]; let i=y[1]; let r=y[2];
        let n=s+i+r;
        dy[0] = -self.beta*s*i/(n+1e-12);
        dy[1] =  self.beta*s*i/(n+1e-12) - self.gamma*i;
        dy[2] =  self.gamma*i;
    }
}
#[wasm_bindgen]
pub fn integrate_sir(beta:f64,gamma:f64,s0:f64,i0:f64,r0:f64,t0:f64,t1:f64,dt:f64)->JsValue{
    let mut p = std::collections::BTreeMap::new();
    p.insert("beta".into(), beta); p.insert("gamma".into(), gamma);
    let series = integrate_model_vec("sir", &p, vec![s0,i0,r0], t0, t1, dt);
    serde_wasm_bindgen::to_value(&series).unwrap()
}


struct Lorenz{ sigma:f64,rho:f64,beta:f64 }
impl System<[f64;3]> for Lorenz {
    fn system(&self,_t:f64,y:&[f64;3],dy:&mut[f64;3]){
        let x=y[0]; let yy=y[1]; let z=y[2];
        dy[0]=self.sigma*(yy-x);
        dy[1]=x*(self.rho - z) - yy;
        dy[2]=x*yy - self.beta*z;
    }
}
#[wasm_bindgen]
pub fn integrate_lorenz63(sigma:f64,rho:f64,beta:f64,x0:f64,y0:f64,z0:f64,t0:f64,t1:f64,dt:f64)->JsValue{
    let mut p = std::collections::BTreeMap::new();
    p.insert("sigma".into(), sigma); p.insert("rho".into(), rho); p.insert("beta".into(), beta);
    let series = integrate_model_vec("lorenz", &p, vec![x0,y0,z0], t0, t1, dt);
    serde_wasm_bindgen::to_value(&series).unwrap()
}


struct Rossler{ a:f64,b:f64,c:f64 }
impl System<[f64;3]> for Rossler {
    fn system(&self,_t:f64,y:&[f64;3],dy:&mut[f64;3]){
        let x=y[0]; let yy=y[1]; let z=y[2];
        dy[0] = -yy - z;
        dy[1] = x + self.a*yy;
        dy[2] = self.b + z*(x - self.c);
    }
}
#[wasm_bindgen]
pub fn integrate_rossler(a:f64,b:f64,c:f64,x0:f64,y0:f64,z0:f64,t0:f64,t1:f64,dt:f64)->JsValue{
    let mut p = std::collections::BTreeMap::new();
    p.insert("a".into(), a); p.insert("b".into(), b); p.insert("c".into(), c);
    let series = integrate_model_vec("rossler", &p, vec![x0,y0,z0], t0, t1, dt);
    serde_wasm_bindgen::to_value(&series).unwrap()
}



// ---------- Simple adaptive RK45 (Cash–Karp) for 3D systems ----------
fn rk45_step_3(f: &dyn Fn(f64, [f64;3]) -> [f64;3], t: f64, y: [f64;3], h: f64) -> ([f64;3], [f64;3]) {
    let a2=0.2;
    let a3=0.3;
    let a4=0.6;
    let a5=1.0;
    let a6=0.875;

    let b21=0.2;

    let b31=3.0/40.0; let b32=9.0/40.0;

    let b41=0.3; let b42=-0.9; let b43=1.2;

    let b51=-11.0/54.0; let b52=2.5; let b53=-70.0/27.0; let b54=35.0/27.0;

    let b61=1631.0/55296.0; let b62=175.0/512.0; let b63=575.0/13824.0; let b64=44275.0/110592.0; let b65=253.0/4096.0;

    let c1=37.0/378.0; let c3=250.0/621.0; let c4=125.0/594.0; let c6=512.0/1771.0;
    let dc1=c1 - 2825.0/27648.0;
    let dc3=c3 - 18575.0/48384.0;
    let dc4=c4 - 13525.0/55296.0;
    let dc5= -277.0/14336.0;
    let dc6=c6 - 0.25;

    let k1 = f(t, y);

    let y2 = [y[0]+h*b21*k1[0], y[1]+h*b21*k1[1], y[2]+h*b21*k1[2]];
    let k2 = f(t + a2*h, y2);

    let y3 = [y[0]+h*(b31*k1[0]+b32*k2[0]), y[1]+h*(b31*k1[1]+b32*k2[1]), y[2]+h*(b31*k1[2]+b32*k2[2])];
    let k3 = f(t + a3*h, y3);

    let y4 = [y[0]+h*(b41*k1[0]+b42*k2[0]+b43*k3[0]), y[1]+h*(b41*k1[1]+b42*k2[1]+b43*k3[1]), y[2]+h*(b41*k1[2]+b42*k2[2]+b43*k3[2])];
    let k4 = f(t + a4*h, y4);

    let y5 = [y[0]+h*(b51*k1[0]+b52*k2[0]+b53*k3[0]+b54*k4[0]), y[1]+h*(b51*k1[1]+b52*k2[1]+b53*k3[1]+b54*k4[1]), y[2]+h*(b51*k1[2]+b52*k2[2]+b53*k3[2]+b54*k4[2])];
    let k5 = f(t + a5*h, y5);

    let y6 = [y[0]+h*(b61*k1[0]+b62*k2[0]+b63*k3[0]+b64*k4[0]+b65*k5[0]),
              y[1]+h*(b61*k1[1]+b62*k2[1]+b63*k3[1]+b64*k4[1]+b65*k5[1]),
              y[2]+h*(b61*k1[2]+b62*k2[2]+b63*k3[2]+b64*k4[2]+b65*k5[2])];
    let k6 = f(t + a6*h, y6);

    let yout = [
        y[0] + h*(c1*k1[0] + c3*k3[0] + c4*k4[0] + c6*k6[0]),
        y[1] + h*(c1*k1[1] + c3*k3[1] + c4*k4[1] + c6*k6[1]),
        y[2] + h*(c1*k1[2] + c3*k3[2] + c4*k4[2] + c6*k6[2]),
    ];

    let err = [
        h*(dc1*k1[0] + dc3*k3[0] + dc4*k4[0] + dc5*k5[0] + dc6*k6[0]),
        h*(dc1*k1[1] + dc3*k3[1] + dc4*k4[1] + dc5*k5[1] + dc6*k6[1]),
        h*(dc1*k1[2] + dc3*k3[2] + dc4*k4[2] + dc5*k5[2] + dc6*k6[2]),
    ];
    (yout, err)
}

fn integrate_adaptive_3(
    f: &dyn Fn(f64, [f64;3]) -> [f64;3],
    t0: f64, t1: f64, mut y: [f64;3],
    dt0: f64, rtol: f64, atol: f64, max_steps: usize
) -> (Vec<f64>, Vec<[f64;3]>) {
    let mut t = t0;
    let mut h = dt0.abs().max(1e-6);
    let mut ts = Vec::new();
    let mut ys = Vec::new();
    ts.push(t); ys.push(y);

    for _ in 0..max_steps {
        if t >= t1 { break; }
        if t + h > t1 { h = t1 - t; }
        let (ynew, err) = rk45_step_3(f, t, y, h);
        let sc0 = atol + rtol * ynew[0].abs().max(y[0].abs());
        let sc1 = atol + rtol * ynew[1].abs().max(y[1].abs());
        let sc2 = atol + rtol * ynew[2].abs().max(y[2].abs());
        let e = ((err[0]/sc0).powi(2) + (err[1]/sc1).powi(2) + (err[2]/sc2).powi(2)).sqrt() / 3.0_f64.sqrt();

        if e <= 1.0 {
            t += h;
            y = ynew;
            ts.push(t);
            ys.push(y);
        }
        let safety = 0.9;
        let pow = 0.2;
        let factor = if e == 0.0 { 5.0 } else { (safety * e.powf(-pow)).clamp(0.2, 5.0) };
        h = (h * factor).clamp(1e-6, 1.0);
    }
    (ts, ys)
}

#[wasm_bindgen]
pub fn integrate_lorenz63_adaptive(sigma:f64,rho:f64,beta:f64,x0:f64,y0:f64,z0:f64,t0:f64,t1:f64,dt0:f64,rtol:f64,atol:f64,max_steps:usize)->JsValue{
    let f = |_:f64, y:[f64;3]| -> [f64;3] {
        let x=y[0]; let yy=y[1]; let z=y[2];
        [sigma*(yy-x), x*(rho-z)-yy, x*yy - beta*z]
    };
    let (ts, ys) = integrate_adaptive_3(&f, t0, t1, [x0,y0,z0], dt0, rtol, atol, max_steps);
    let series = Series{ t: ts, y: flatten(&ys), dim:3 };
    serde_wasm_bindgen::to_value(&series).unwrap()
}

#[wasm_bindgen]
pub fn integrate_rossler_adaptive(a:f64,b:f64,c:f64,x0:f64,y0:f64,z0:f64,t0:f64,t1:f64,dt0:f64,rtol:f64,atol:f64,max_steps:usize)->JsValue{
    let f = |_:f64, y:[f64;3]| -> [f64;3] {
        let x=y[0]; let yy=y[1]; let z=y[2];
        [-yy - z, x + a*yy, b + z*(x - c)]
    };
    let (ts, ys) = integrate_adaptive_3(&f, t0, t1, [x0,y0,z0], dt0, rtol, atol, max_steps);
    let series = Series{ t: ts, y: flatten(&ys), dim:3 };
    serde_wasm_bindgen::to_value(&series).unwrap()
}

// ---------- Extra advanced systems (fixed-step) ----------
struct Chua{ alpha:f64,beta:f64,m0:f64,m1:f64 }
impl System<[f64;3]> for Chua {
    fn system(&self,_t:f64,y:&[f64;3],dy:&mut[f64;3]){
        let x=y[0]; let yy=y[1]; let z=y[2];
        let fx = self.m1*x + 0.5*(self.m0-self.m1)*((x+1.0).abs() - (x-1.0).abs());
        dy[0]=self.alpha*(yy - x - fx);
        dy[1]=x - yy + z;
        dy[2]=-self.beta*yy;
    }
}
#[wasm_bindgen]
pub fn integrate_chua(alpha:f64,beta:f64,m0:f64,m1:f64,x0:f64,y0:f64,z0:f64,t0:f64,t1:f64,dt:f64)->JsValue{
    let mut p = std::collections::BTreeMap::new();
    p.insert("alpha".into(), alpha); p.insert("beta".into(), beta);
    p.insert("m0".into(), m0); p.insert("m1".into(), m1);
    let series = integrate_model_vec("chua", &p, vec![x0,y0,z0], t0, t1, dt);
    serde_wasm_bindgen::to_value(&series).unwrap()
}


struct Thomas{ a:f64 }
impl System<[f64;3]> for Thomas {
    fn system(&self,_t:f64,y:&[f64;3],dy:&mut[f64;3]){
        let x=y[0]; let yy=y[1]; let z=y[2];
        dy[0]=yy.sin() - self.a*x;
        dy[1]=z.sin()  - self.a*yy;
        dy[2]=x.sin()  - self.a*z;
    }
}
#[wasm_bindgen]
pub fn integrate_thomas(a:f64,x0:f64,y0:f64,z0:f64,t0:f64,t1:f64,dt:f64)->JsValue{
    let mut p = std::collections::BTreeMap::new();
    p.insert("a".into(), a);
    let series = integrate_model_vec("thomas", &p, vec![x0,y0,z0], t0, t1, dt);
    serde_wasm_bindgen::to_value(&series).unwrap()
}

