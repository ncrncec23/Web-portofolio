class Persp {
    constructor(Platno, Xmin, Xmax, ry, d) {
        this.platno = Platno;
        this.xmin = Xmin;
        this.xmax = Xmax;
        const width = this.platno.width;
        const height = this.platno.height;
        this.d = d;
        this.sx = width / (this.xmax - this.xmin);
        this.sy = -this.sx;
        this.px = -this.sx * this.xmin;
        this.py = (1 - ry) * height;
        this.ymax = this.py / this.sx;
        this.ymin = this.ymax - height / this.sx;
        this.M = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
        this.stogM = [];
        this.g = this.platno.getContext("2d");
    }

    pravokutnik(a, b, boja) {
        this.koristiBoju(boja);
        this.postaviNa(0, 0);
        this.linijaDo(a, 0);
        this.linijaDo(a, b);
        this.linijaDo(a - a, b);
        this.linijaDo(a - a, b - b);
        this.povuciLiniju();
    }


    spremi() {
        this.stogM.push(this.M);
    }

    vrati() {
        this.M = this.stogM.pop();
    }

    //Metoda moveTo();
    postaviNa(x, y, z = 0) {
        this.g.beginPath();
        let [xc, yc] = this.TK(x, y, z);
        this.g.moveTo(xc, yc);
    }

    //Metoda lineTo();
    linijaDo(x, y, z = 0) {
        let [xc, yc] = this.TK(x, y, z);
        this.g.lineTo(xc, yc);
    }

    //Metoda boje
    koristiBoju(c) {
        this.g.strokeStyle = c;
    }

    //Metoda za stroke();
    povuciLiniju() {
        this.g.stroke();
    }

    sinus() {
        this.koristiBoju("red");
        for (let n = 0; n <= 12; n++) {
            this.postaviNa(0, 0);
            for (let i = 0; i <= (2 * Math.PI); i += 0.01) {
                this.linijaDo(i, Math.sin(i));
            }
            this.povuciLiniju();
            this.rotiraj(n * 30);
        }
    }

    kruznica(R, nseg, boja) {
        const dphi = 2 * Math.PI / nseg;
        this.postaviNa(R, 0);
        this.koristiBoju(boja);
        for (let i = 1; i <= nseg; i++) {
            let x = R * Math.cos(i * dphi);
            let y = R * Math.sin(i * dphi);
            this.linijaDo(x, y);
        }
        this.povuciLiniju();
    }

    elipsa(a, b, nseg, boja) {
        const dphi = 2 * Math.PI / nseg;
        this.postaviNa(a, 0);
        this.koristiBoju(boja);
        for (let i = 1; i <= nseg; i++) {
            let x = a * Math.cos(i * dphi);
            let y = b * Math.sin(i * dphi);
            this.linijaDo(x, y);
        }
        this.povuciLiniju();
    }

    puta(mat) {
        let rezultat = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++)
                for (let k = 0; k < 4; k++)
                    rezultat[i][j] += this.M[i][k] * mat[k][j];
        }
        this.M = rezultat;
    }

    skaliraj(sx, sy, sz = 0) {
        this.puta([[sx, 0, 0, 0], [0, sy, 0, 0], [0, 0, sz, 0], [0, 0, 0, 1]]);
    }

    pomakni(px, py, pz = 0) {
        this.puta([[1, 0, 0, px], [0, 1, 0, py], [0, 0, 1, pz], [0, 0, 0, 1]]);
    }

    rotirajX(kut) {
        let phi = kut * Math.PI / 180;
        this.puta(
            [[1, 0, 0, 0],
            [0, Math.cos(phi), -Math.sin(phi), 0],
            [0, Math.sin(phi), Math.cos(phi), 0],
            [0, 0, 0, 1]]);
    } // rotirajX

    rotirajY(kut) {
        let phi = kut * Math.PI / 180;
        this.puta(
            [[Math.cos(phi), 0, Math.sin(phi), 0],
            [0, 1, 0, 0],
            [-Math.sin(phi), 0, Math.cos(phi), 0],
            [0, 0, 0, 1]]);
    } // rotirajY

    rotirajZ(kut) {
        let phi = kut * Math.PI / 180;
        this.puta(
            [[Math.cos(phi), -Math.sin(phi), 0, 0],
            [Math.sin(phi), Math.cos(phi), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]]);
    } // rotirajZ

    rotiraj(kut) {
        let phi = (kut * Math.PI) / 180;
        this.puta([[Math.cos(phi), -Math.sin(phi), 0], [Math.sin(phi), Math.cos(phi), 0], [0, 0, 1]]);
    }

    identitet() {
        this.M = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    }

    TK(x, y, z) {
        // transformacija x i y u globalne koordinate:
        // množenje s matricom transformacije this.M
        let xg = this.M[0][0] * x + this.M[0][1] * y + this.M[0][2] * z + this.M[0][3];
        let yg = this.M[1][0] * x + this.M[1][1] * y + this.M[1][2] * z + this.M[1][3];
        let zg = this.M[2][0] * x + this.M[2][1] * y + this.M[2][2] * z + this.M[2][3];

        let xt = -this.d * xg / zg;
        let yt = -this.d * yg / zg;
        // xg i yg preračunati u zaslonske koordinate xc i yc
        let xc = this.sx * xt + this.px;
        let yc = this.sy * yt + this.py;
        return [xc, yc];
    }

    pravac(k, l, boja) {
        this.koristiBoju(boja);
        let yl = k * this.xmin + l;
        let yd = k * this.xmax + l;
        this.postaviNa(this.xmin, yl);
        this.linijaDo(this.xmax, yd);
        this.povuciLiniju();
    }

    rotirajOko(u, v, kut) {
        this.pomakni(u, v);
        this.rotiraj(kut);
        this.pomakni(-u, -v);
    }

    iscrtajOsi() {
        this.koristiBoju("black");
        this.postaviNa(this.xmin, 0);
        this.linijaDo(this.xmax, 0);
        this.povuciLiniju();
        this.postaviNa(this.xmin, 0);
        this.linijaDo(this.xmax, 0);
        this.postaviNa(0, this.ymax);
        this.linijaDo(0, this.ymin);
        this.povuciLiniju();
        for (let i = 0; i <= this.xmax; i++) {
            this.postaviNa(i, 1 / 10);
            this.linijaDo(i, -(1 / 10));
            this.povuciLiniju();
        }
        for (let i = 0; i >= this.xmin; i--) {
            this.postaviNa(i, 1 / 10);
            this.linijaDo(i, -(1 / 10));
            this.povuciLiniju();
        }
        for (let i = 0; i <= this.ymax; i++) {
            this.postaviNa(1 / 10, i);
            this.linijaDo(-(1 / 10), i);
            this.povuciLiniju();
        }
        for (let i = 0; i >= this.ymin; i--) {
            this.postaviNa(1 / 10, i);
            this.linijaDo(-(1 / 10), i);
            this.povuciLiniju();
        }

    }
    osi3D(l) {
        this.koristiBoju("red");
        this.postaviNa(0, 0);
        this.linijaDo(l, 0);
        this.povuciLiniju();
        this.koristiBoju("green");
        this.postaviNa(0, 0);
        this.linijaDo(0, l);
        this.povuciLiniju();
        this.koristiBoju("blue");
        this.postaviNa(0, 0, 0);
        this.linijaDo(0, 0, l);
        this.povuciLiniju();
    }

    kocka(a, boja) {
        this.koristiBoju(boja);
        this.postaviNa(0, 0);
        this.linijaDo(0, a);
        this.linijaDo(a, a);
        this.linijaDo(a, a - a);
        this.linijaDo(a - a, a - a);
        this.povuciLiniju();
        this.postaviNa(0, 0, a);
        this.linijaDo(0, a, a);
        this.linijaDo(a, a, a);
        this.linijaDo(a, a - a, a);
        this.linijaDo(a - a, a - a, a);
        this.povuciLiniju();
        this.postaviNa(0, 0, a);
        this.linijaDo(0, 0, a - a);
        this.povuciLiniju();
        this.postaviNa(0, a, a);
        this.linijaDo(0, a, a - a);
        this.povuciLiniju();
        this.postaviNa(a, a, a);
        this.linijaDo(a, a, a - a);
        this.povuciLiniju();
        this.postaviNa(a, a - a, a);
        this.linijaDo(a, a - a, a - a);
        this.povuciLiniju();
    }

    postaviKameru(x0, y0, z0, x1, y1, z1, Vx, Vy, Vz) {
        function normiraj(v) {
            let L = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            return [v[0] / L, v[1] / L, v[2] / L];
        } // normiraj

        function vektorskiProdukt(ax, ay, az, bx, by, bz) {
            return [ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx];
        } // vektorskiProdukt

        let N = [x0 - x1, y0 - y1, z0 - z1];
        let n = normiraj(N);
        let U = vektorskiProdukt(Vx, Vy, Vz, n[0], n[1], n[2]);
        let u = normiraj(U);
        let v = vektorskiProdukt(n[0], n[1], n[2], u[0], u[1], u[2]);

        this.puta(
            [[u[0], u[1], u[2], -u[0] * x0 - u[1] * y0 - u[2] * z0],
            [v[0], v[1], v[2], -v[0] * x0 - v[1] * y0 - v[2] * z0],
            [n[0], n[1], n[2], -n[0] * x0 - n[1] * y0 - n[2] * z0],
            [0, 0, 0, 1]]
        )
    }

    zelenamreza() {
        this.pomakni(-10, -10)
        this.koristiBoju('00000052');
        for (let i = 0; i <= 20; i++) {
            this.postaviNa(0, i);
            this.linijaDo(20, i);
            this.povuciLiniju();
        }
        for (let i = 0; i <= 20; i++) {
            this.postaviNa(i, 0);
            this.linijaDo(i, 20);
            this.povuciLiniju();
        }
        this.pomakni(10, 10)
    }

    stozac(R, nseg, h, boja) {
        this.kruznica(R, nseg, boja)
        const dphi = 2 * Math.PI / nseg;
        this.postaviNa(R, 0);
        this.koristiBoju(boja);
        for (let i = 1; i <= nseg; i++) {
            let x = R * Math.cos(i * dphi);
            let y = R * Math.sin(i * dphi);
            this.linijaDo(x, y);
            this.linijaDo(0, 0, h);
        }
        this.povuciLiniju();
    }

    valjak(R, nseg, h, boja) {
        this.kruznica(R, nseg, boja)
        this.pomakni(0, 0, h)
        this.kruznica(R, nseg, boja)
        this.pomakni(0, 0, -h)
        const dphi = 2 * Math.PI / nseg;
        this.postaviNa(R, 0);
        this.koristiBoju(boja);
        for (let i = 1; i <= nseg; i++) {
            let x = R * Math.cos(i * dphi);
            let y = R * Math.sin(i * dphi);
            this.postaviNa(x, y);
            this.linijaDo(x, y, h);
            this.povuciLiniju();
        }
    }

    kugla(R, m, n, boja) {
        let nseg = 50;
        const dphi = 2 * Math.PI / nseg;
        this.koristiBoju(boja);
        //crtanje paralela
        for (let teta = 0; teta <= Math.PI; teta += Math.PI / (n + 1)) {
            //let p = R * Math.sin(teta);
            this.postaviNa(R * Math.sin(teta) * Math.cos(dphi), R * Math.sin(teta) * Math.sin(dphi), R * Math.cos(teta));
            for (let i = 0; i <= nseg; i++) {
                let x = R * Math.sin(teta) * Math.cos(i * dphi);
                //let x = p * Math.cos(i * dphi);
                let y = R * Math.sin(teta) * Math.sin(i * dphi);
                //let y = p * Math.sin(i * dphi);
                let z = R * Math.cos(teta)
                this.linijaDo(x, y, z);
            }
            this.povuciLiniju();
        }

        //crtanje meridijana
        for (let phi = 0; phi <= 2 * Math.PI; phi += 2 * Math.PI / m) {
            this.postaviNa(0, 0, R);
            for (let teta = 0; teta <= Math.PI; teta += 0.01) {
                let x = R * Math.sin(teta) * Math.cos(phi);
                let y = R * Math.sin(teta) * Math.sin(phi);
                let z = R * Math.cos(teta);
                this.linijaDo(x, y, z);
            }
            this.povuciLiniju();
        }

        //ISCRTAVANJE SAMO JEDNOG MERIDIJANA
        /*this.postaviNa(0, 0, R);
        for (let kut = 0; kut <= Math.PI; kut += 0.01) {
            let x = R * Math.sin(kut);
            let z = R * Math.cos(kut);
            this.linijaDo(x, 0, z);
        }
        this.povuciLiniju();
        */
    }

    krnjistozac(R, R2, nseg, h, boja) {
        this.kruznica(R, nseg, boja)
        const dphi = 2 * Math.PI / nseg;
        this.postaviNa(R, 0);
        this.koristiBoju(boja);
        for (let i = 1; i <= nseg; i++) {
            let x = R * Math.cos(i * dphi);
            let y = R * Math.sin(i * dphi);
            let x2 = R2 * Math.cos(i * dphi);
            let y2 = R2 * Math.sin(i * dphi);
            this.postaviNa(x, y);
            this.linijaDo(x2, y2, h);
            this.povuciLiniju();
        }
        this.pomakni(0, 0, h);
        this.kruznica(R2, nseg, boja);
    }

    polukugla(R, m, n, boja) {
        let nseg = 50;
        const dphi = 2 * Math.PI / nseg;
        this.koristiBoju(boja);
        //crtanje paralela
        for (let teta = 0; teta <= Math.PI / 2; teta += Math.PI / (n + 1)) {
            //let p = R * Math.sin(teta);
            this.postaviNa(R * Math.sin(teta) * Math.cos(dphi), R * Math.sin(teta) * Math.sin(dphi), R * Math.cos(teta));
            for (let i = 0; i <= nseg; i++) {
                let x = R * Math.sin(teta) * Math.cos(i * dphi);
                //let x = p * Math.cos(i * dphi);
                let y = R * Math.sin(teta) * Math.sin(i * dphi);
                //let y = p * Math.sin(i * dphi);
                let z = R * Math.cos(teta)
                this.linijaDo(x, y, z);
            }
            this.povuciLiniju();
        }
        //crtanje meridijana
        for (let phi = 0; phi <= 2 * Math.PI; phi += 2 * Math.PI / m) {
            this.postaviNa(0, 0, R);
            for (let teta = 0; teta <= Math.PI / 2; teta += 0.01) {
                let x = R * Math.sin(teta) * Math.cos(phi);
                let y = R * Math.sin(teta) * Math.sin(phi);
                let z = R * Math.cos(teta);
                this.linijaDo(x, y, z);
            }
            this.povuciLiniju();
        }
    }
}

function mobileLanchange() {
    const gumbOne = document.getElementById("nemaMargin1");
    const gumbThree = document.getElementById("nemaMargin3");
    const gumbFive = document.getElementById("nemaMargin5");
    gumbOne.addEventListener("click", () => changeLanguage("hr"));
    gumbThree.addEventListener("click", () => changeLanguage("en"));
    gumbFive.addEventListener("click", () => changeLanguage("de"));
}

function mobileAccesibility() {
    const gumb = document.querySelectorAll(".unut");
    gumb.forEach((el) => {
        el.onclick = () => {
            if (el.classList.contains("unutarnji_klik")) {
                el.classList.remove("unutarnji_klik");
                el.classList.add("unatrag_unutarnji_klik");
            } else {
                el.classList.remove("unatrag_unutarnji_klik");
                el.classList.add("unutarnji_klik");
            }
        }
    });
}

function theme_mobile() {
    const svege1 = document.getElementById("ESVEGE1");
    const svege2 = document.getElementById("ESVEGE2");
    const path1 = document.getElementById("PATH1");
    const path2 = document.getElementById("PATH2");
    svege1.onclick = () => {
        path1.style.fill = "#2E1A47";
        path2.style.fill = "#8C898E";
    }
    svege2.onclick = () => {
        path1.style.fill = "#8C898E";
        path2.style.fill = "#F7CC74";
    }
}

const button = document.getElementById("setts");
const dropdown = document.querySelector(".dropdownSettings");
const dropdownP = document.querySelectorAll(".settingsText");
const but = document.querySelectorAll(".but");
const svg = document.getElementById("svgSetts");
const seperator = document.querySelectorAll(".seperator");

const language = document.getElementById("Language");
const accesibility = document.getElementById("Accessibility");
const theme = document.getElementById("Theme");

const postavke = document.querySelector(".settings");
const mq = window.matchMedia("(max-width: 768px)");

function provjeriMediaQuery(e) {
    mobileLanchange();
    mobileAccesibility();
    theme_mobile();
    const settingsButon = document.querySelector(".settingsBut");
    const latoBold = document.querySelectorAll(".lato-bold");
    if (e.matches) {
        settingsButon.style.display = "none";
        postavke.innerHTML += `
            <div class="HamburgerContainer">
                <div class="firstLineHam"></div>
                <div class="secondLineHam"></div>
                <div class="thirdLineHam"></div>
            </div>`
        latoBold.forEach((e) => {
            e.style.display = "none";
        });
        let open = false;
        const hamburger = document.querySelector(".HamburgerContainer");
        const hamburger2 = document.querySelector(".HamburgerContainer2");
        const prvaLinija2 = document.querySelector(".firstLineHam2");
        const trecaLinija2 = document.querySelector(".thirdLineHam2");
        const drugaLinija2 = document.querySelector(".secondLineHam2");
        const prvaLinija = document.querySelector(".firstLineHam");
        const trecaLinija = document.querySelector(".thirdLineHam");
        const drugaLinija = document.querySelector(".secondLineHam");
        const sideBar = document.querySelector(".sidebar");
        const mobileNav = document.querySelector(".mobileNavigation");
        const mobilnaNav = document.querySelector(".mobilnaNav");
        const LAN = document.getElementById("LAN");
        const acc = document.getElementById("acc");
        const swich = document.getElementById("switch");
        const tm = document.getElementById("tm");
        const drklght = document.querySelector(".dark_light");
        const black = document.querySelector(".black");
        const JeziciPromjena = document.querySelector(".jeziciPromjena");
        hamburger.onclick = () => {
            if (!open) {
                prvaLinija.classList.remove("klikPrvaReverse");
                trecaLinija.classList.remove("klikTrecaReverse");
                drugaLinija.classList.remove("transparityReverse");
                prvaLinija.classList.add("klikPrva");
                trecaLinija.classList.add("klikTreca");
                drugaLinija.classList.add("transparity");
                prvaLinija2.classList.remove("klikPrvaReverse2");
                trecaLinija2.classList.remove("klikTrecaReverse2");
                drugaLinija2.classList.remove("transparityReverse2");
                prvaLinija2.classList.add("klikPrva2");
                trecaLinija2.classList.add("klikTreca2");
                drugaLinija2.classList.add("transparity2");
                mobileNav.classList.remove("flexingReverse");
                mobilnaNav.classList.remove("flexingReverse");
                mobileNav.style.display = "flex";
                mobileNav.classList.add("flexing");
                mobilnaNav.classList.add("flexing");
                black.classList.add("flexing");
                sideBar.classList.remove("sidebarClose");
                sideBar.classList.add("sidebarOpen");
                if (sideBar.classList.contains("sidebarOpen")) {
                    document.body.style.overflow = "hidden";
                }
            } else {
                sideBar.classList.remove("sidebarOpen");
                sideBar.classList.add("sidebarClose");
                prvaLinija.classList.remove("klikPrva");
                trecaLinija.classList.remove("klikTreca");
                drugaLinija.classList.remove("transparity");
                prvaLinija.classList.add("klikPrvaReverse");
                trecaLinija.classList.add("klikTrecaReverse");
                drugaLinija.classList.add("transparityReverse");
            }
            open = !open
        }

        hamburger2.onclick = () => {
            black.classList.remove("flexing");
            mobileNav.classList.remove("flexing");
            mobilnaNav.classList.remove("flexing");
            mobileNav.classList.add("flexingReverse");
            mobilnaNav.classList.add("flexingReverse");
            black.classList.add("flexingReverse");
            dropdown.classList.remove("show");
            dropdown.classList.add("unshow");
            prvaLinija.classList.remove("klikPrva");
            trecaLinija.classList.remove("klikTreca");
            drugaLinija.classList.remove("transparity");
            prvaLinija.classList.add("klikPrvaReverse");
            trecaLinija.classList.add("klikTrecaReverse");
            drugaLinija.classList.add("transparityReverse");
            prvaLinija2.classList.remove("klikPrva2");
            trecaLinija2.classList.remove("klikTreca2");
            drugaLinija2.classList.remove("transparity2");
            prvaLinija2.classList.add("klikPrvaReverse2");
            trecaLinija2.classList.add("klikTrecaReverse2");
            drugaLinija2.classList.add("transparityReverse2");
            sideBar.classList.remove("sidebarOpen");
            sideBar.classList.add("sidebarClose");
            if (sideBar.classList.contains("sidebarClose")) {
                document.body.style.overflow = "scroll";
            }
            open = false;
        }

        const transformers = document.querySelectorAll(".transformers");

        LAN.onclick = () => {
            if (JeziciPromjena.classList.contains("otkriveno")) {
                JeziciPromjena.classList.remove("otkriveno");
                JeziciPromjena.classList.add("sakriveno");
                transformers[0].classList.add("k");
                transformers[0].classList.remove("kr");
            } else {
                JeziciPromjena.classList.remove("sakriveno");
                JeziciPromjena.classList.add("otkriveno");
                transformers[0].classList.add("kr");
                transformers[0].classList.remove("k");
            }
        }

        acc.onclick = () => {
            if (swich.classList.contains("opn")) {
                swich.classList.remove("opn");
                swich.classList.add("cls");
                transformers[1].classList.remove("kr");
                transformers[1].classList.add("k");
            } else {
                swich.classList.remove("cls");
                swich.classList.add("opn");
                transformers[1].classList.remove("k");
                transformers[1].classList.add("kr");
            }
        }

        tm.onclick = () => {
            if (drklght.classList.contains("o")) {
                drklght.classList.remove("o");
                drklght.classList.add("c");
                transformers[2].classList.remove("kr");
                transformers[2].classList.add("k");
            } else {
                drklght.classList.remove("c");
                drklght.classList.add("o");
                transformers[2].classList.remove("k");
                transformers[2].classList.add("kr");
            }
        }

    } else {
        settingsButon.style.display = "inline-block";
        latoBold.forEach((e) => {
            e.style.display = "inline-block";
        });
        const existingHamburger = document.querySelector(".HamburgerContainer");
        if (existingHamburger) {
            existingHamburger.remove();
        }
    }

}

provjeriMediaQuery(mq);

mq.addEventListener("change", provjeriMediaQuery);

const c = document.getElementById("canvas");
const p = document.getElementById("canvas1");

let crtanje = new Persp(c, -16, 16, 1 / 2, 18);
let crt = new Persp(p, -16, 16, 1 / 2, 18);

const ctx = c.getContext("2d");
const ctxx = p.getContext("2d");

const w = c.width;
const h = c.height;
const pw = p.width;
const ph = p.height;

let kut = 0;
let kutt = 0;

window.onload = () => {
    paint();
    iscrtaj();
}

function paint() {
    ctxx.fillStyle = "#213555";
    ctxx.fillRect(0, 0, pw, ph);
    crt.postaviKameru(20, 20, 5, 0, 0, 0, 0, 0, 1);
    crt.rotirajX(kut)
    crt.rotirajY(-50)
    crt.kugla(15, 40, 25, "#00000060");
    crt.identitet();
    kut++;
    requestAnimationFrame(paint);
}


function iscrtaj() {
    ctx.fillStyle = "#213555";
    ctx.fillRect(0, 0, w, h);
    crtanje.postaviKameru(14, 14, 19, 0, 0, 0, 0, 0, 1);
    crtanje.zelenamreza();
    crtanje.stozac(5, 25, 10, "#00000052");
    crtanje.rotirajZ(kutt)
    crtanje.pomakni(0, 0, 8.5)
    crtanje.valjak(1, 20, 2, "#00000052")
    crtanje.spremi();
    crtanje.spremi();
    crtanje.pomakni(1, 0, 1)
    crtanje.rotirajY(90);
    crtanje.valjak(0.5, 10, 6, "#00000052")
    crtanje.pomakni(0, 0, 9)
    crtanje.rotirajX(-90);
    crtanje.polukugla(3, 32, 15, "#00000052");
    crtanje.vrati();
    crtanje.pomakni(-0.45, -0.85, 1)
    crtanje.rotirajZ(-30);
    crtanje.rotirajX(90);
    crtanje.valjak(0.5, 10, 6, "#00000052")
    crtanje.pomakni(0, 0, 9)
    crtanje.rotirajY(90);
    crtanje.polukugla(3, 32, 15, "#00000052");
    crtanje.vrati();
    crtanje.pomakni(-0.5, 0.82, 1)
    crtanje.rotirajZ(30);
    crtanje.rotirajX(-90);
    crtanje.valjak(0.5, 10, 6, "#00000052")
    crtanje.pomakni(0, 0, 9)
    crtanje.rotirajY(-90);
    crtanje.polukugla(3, 32, 15, "#00000052");
    crtanje.identitet();
    kutt++;
    requestAnimationFrame(iscrtaj)
}

const roditelj = document.querySelector(".content");
const dijete = document.createElement("div");
dijete.classList.add("dijete");

let placeHolder = document.createElement("p");
placeHolder.textContent = "";
placeHolder.style.margin = "0";
dijete.appendChild(placeHolder)
roditelj.appendChild(dijete);

function language_change() {

    if (dijete.hasChildNodes() && roditelj.hasChildNodes()) {
        let hr = document.createElement("a");
        hr.setAttribute("href", "index.html")
        hr.textContent = "HR"
        let sep = document.createElement("div");
        sep.textContent = "/";
        let en = document.createElement("a");
        en.setAttribute("href", "index.html");
        en.textContent = "EN"
        let sep2 = document.createElement("div");
        sep2.textContent = "/";
        let de = document.createElement("a");
        de.setAttribute("href", "index.html");
        de.textContent = "DE";
        dijete.replaceChildren(hr, sep, en, sep2, de);
        roditelj.replaceChildren(dijete);
        hr.addEventListener("click", () => changeLanguage("hr"));
        en.addEventListener("click", () => changeLanguage("en"));
        de.addEventListener("click", () => changeLanguage("de"));
    }
}

function lan_change() {
    let savedLang = localStorage.getItem("language");
    if (savedLang) {
        changeLanguage(savedLang);
    }
}

function changeLanguage(lang) {
    const home = document.getElementById("home");
    const aboutme = document.getElementById("aboutme");
    const projects = document.getElementById("pj");
    const helloEveryone = document.getElementById("helloEveryone");
    const mojeIme = document.getElementById("mojeImeje");
    const lm = document.getElementById("lm");
    const edu = document.getElementById("edu");
    const opisedu = document.getElementById("opisEdu");
    const procitaj = document.getElementById("procitajvise");
    const wrk = document.getElementById("wrk");
    const wrkopis = document.getElementById("wrkopis");
    const rdm = document.getElementById("rdm");
    const vještine = document.getElementById("vještine");
    const vještineopis = document.getElementById("vještineopis");
    const rdmvještine = document.getElementById("rdmvještine");
    const abtm = document.getElementById("abtm");
    const opss = document.getElementById("opss");
    const svaPrava = document.getElementById("svaPrava");
    const DevNiko = document.getElementById("DevNiko");
    const fotA = document.getElementById("fotA");
    const fotP = document.getElementById("fotP");
    const fotH = document.getElementById("fotH");
    const javise = document.getElementById("javise");
    const ž = document.getElementById("ž");
    const lj = document.getElementById("lj");
    const nj = document.getElementById("nj");
    document.documentElement.lang = lang;
    localStorage.setItem("language", lang);
    if (lang == "hr") {
        ž.textContent = "Jezik"
        lj.textContent = "Pristupačnost"
        nj.textContent = "Tema"
        home.textContent = "Naslovnica"
        aboutme.textContent = "O meni"
        projects.textContent = "Projekti"
        helloEveryone.textContent = "Pozdrav svima,"
        mojeIme.innerHTML =
            `Moje ime je <span class="lato-bold-p">Niko Crnčec</span>, student sam Fakulteta organizacije i informatike u Varaždinu. 
        Moji hobiji su <span class="lato-bold-p">programiranje</span> i <span class="lato-bold-p">sport</span>. 
        Volim razvijati web, stolne i mobilne aplikacije.`
        lm.textContent = "Saznaj više"
        edu.textContent = "Obrazovanje"
        opisedu.textContent = "Završio sam osnovnu školu Sračinec, zatim Elektrostrojarsku školu Varaždin. Sada studiram poslovnu informatiku na FOI-u";
        procitaj.textContent = "Pročitaj više";
        wrk.textContent = "Posao";
        wrkopis.textContent = "Trenutno nemam radnog iskustva u informacijskim tehnologijama, međutim radio sam kao blagajnik u Lidlu";
        rdm.textContent = "Pročitaj više";
        vještine.textContent = "Vještine";
        vještineopis.textContent = `Imam iskustva s HTML-om, CSS-om i JavaScriptom za razvoj web aplikacija. Također, vješt sam u radu s bazama podataka, kao i u programiranju u C++ i C#`
        rdmvještine.textContent = "Pročitaj više";
        abtm.textContent = "O MENI"
        opss.textContent = `Rođen sam u Varaždinu, Hrvatska, a cijeli svoj život živim u Svibovcu Podravskom. 
        Mojim strastima pripadaju programiranje i sport, a uživam u dobroj hrani i čaši crnog vina. U budućnosti se vidim u uspješnoj IT organizaciji 
        koja razvija inovativna rješenja za velika poduzeća. Smatram se skromnom i poniznom osobom, 
        spremnom na nove izazove, te čekam priliku koja će mi omogućiti da unaprijedim svoj život i ostvarim svoje ciljeve.`
        svaPrava.innerHTML = "&copy; 2025 Niko Crnčec. Sva prava pridržana.";
        DevNiko.textContent = "Razvijeno od strane Niko Crnčec.";
        fotA.textContent = "O meni";
        fotP.textContent = "Projekti";
        fotH.textContent = "Naslovnica";
        javise.textContent = "Kontaktiraj me";
    } else if (lang == "en") {
        //podeseno default
    } else {

    }
}

document.addEventListener("DOMContentLoaded", lan_change);

function theme_change() {
    if (dijete.hasChildNodes() && roditelj.hasChildNodes()) {
        dijete.innerHTML = "";
        dijete.innerHTML += `
        <svg xmlns="http://www.w3.org/2000/svg" width="22" height="22" viewBox="0 0 22 22" fill="none" id="svg" class="opala">
        <path id="pat" d="M11.0778 21.6665C8.11513 21.6665 5.59647 20.6292 3.5218 18.5545C1.44713 16.4799 0.410244 13.9616 0.411133 10.9999C0.411133 8.80431 1.03913 6.78875 2.29513 4.95319C3.55202 3.11586 5.33424 1.79453 7.6418 0.989194C7.9138 0.894083 8.15202 0.869194 8.35647 0.914527C8.56091 0.959861 8.73024 1.04831 8.86447 1.17986C8.99869 1.31142 9.08402 1.48031 9.12047 1.68653C9.15602 1.89364 9.12802 2.11097 9.03647 2.33853C8.86402 2.76253 8.73869 3.19408 8.66047 3.63319C8.58224 4.07231 8.54358 4.52786 8.54447 4.99986C8.54447 7.37586 9.37247 9.39231 11.0285 11.0492C12.6854 12.7052 14.7018 13.5332 17.0778 13.5332C17.6982 13.5332 18.2658 13.4674 18.7805 13.3359C19.296 13.2043 19.7351 13.0976 20.0978 13.0159C20.2916 12.9803 20.4689 12.9848 20.6298 13.0292C20.7907 13.0736 20.9205 13.1519 21.0191 13.2639C21.1214 13.375 21.1907 13.5119 21.2271 13.6745C21.2636 13.8372 21.2427 14.0225 21.1645 14.2305C20.5289 16.4083 19.2831 18.1919 17.4271 19.5812C15.5711 20.9705 13.4547 21.6656 11.0778 21.6665Z" fill="#8C898E"/>
        </svg>`
        dijete.innerHTML += `<div class="opala">/</div>`;
        dijete.innerHTML += `
        <svg xmlns="http://www.w3.org/2000/svg" width="28" height="28" viewBox="0 0 28 28" fill="none" id="svgL" class="opala">
        <path id="patL" d="M13.9998 23.3332C14.3264 23.3332 14.6416 23.4531 14.8857 23.6701C15.1297 23.8871 15.2856 24.1862 15.3238 24.5105L15.3332 
        24.6665V25.9998C15.3328 26.3397 15.2027 26.6665 14.9694 26.9137C14.7361 27.1608 14.4172 27.3095 14.078 27.3294C13.7387 27.3493 13.4047 27.2389 
        13.144 27.0208C12.8834 26.8027 12.716 26.4933 12.6758 26.1558L12.6665 25.9998V24.6665C12.6665 24.3129 12.807 23.9737 13.057 23.7237C13.3071 
        23.4736 13.6462 23.3332 13.9998 23.3332ZM22.4172 20.5465L22.5425 20.6572L23.4758 21.5905C23.715 21.8305 23.8538 22.1524 23.8641 22.491C23.8745 
        22.8296 23.7555 23.1595 23.5315 23.4135C23.3074 23.6676 22.9951 23.8269 22.6579 23.859C22.3206 23.8912 21.9838 23.7937 21.7158 23.5865L21.5905 
        23.4758L20.6572 22.5425C20.4271 22.3128 20.2888 22.007 20.2684 21.6825C20.248 21.358 20.3469 21.0373 20.5464 20.7806C20.7459 20.5239 21.0324 
        20.349 21.3518 20.2887C21.6713 20.2284 22.0018 20.2869 22.2812 20.4532L22.4172 20.5465ZM7.34249 20.6572C7.57206 20.8868 7.70997 21.1922 7.73034 21.5163C7.75072 21.8403 7.65216 22.1606 7.45316 22.4172L7.34249 22.5425L6.40916 23.4758C6.16922 23.715 5.84725 23.8538 5.50865 23.8641C5.17005 23.8745 4.84021 23.7556 4.58612 23.5315C4.33203 23.3075 4.17275 22.9951 4.14063 22.6579C4.1085 22.3206 4.20594 21.9838 4.41316 21.7158L4.52383 21.5905L5.45716 20.6572C5.70719 20.4072 6.04627 20.2668 6.39983 20.2668C6.75338 20.2668 7.09246 20.4072 7.34249 20.6572ZM3.33316 12.6665C3.673 12.6669 3.99987 12.797 4.24698 13.0303C4.4941 13.2636 4.6428 13.5824 4.66272 13.9217C4.68264 14.261 4.57226 14.595 4.35414 14.8556C4.13602 15.1162 3.82662 15.2837 3.48916 15.3238L3.33316 15.3332H1.99983C1.65999 15.3328 1.33312 15.2027 1.086 14.9694C0.838886 14.7361 0.690178 14.4172 0.670262 14.078C0.650346 13.7387 0.760724 13.4047 0.978845 13.1441C1.19697 12.8835 1.50636 12.716 1.84382 12.6758L1.99983 12.6665H3.33316ZM25.9998 12.6665C26.3397 12.6669 26.6665 12.797 26.9136 13.0303C27.1608 13.2636 27.3095 13.5824 27.3294 13.9217C27.3493 14.261 27.2389 14.595 27.0208 14.8556C26.8027 15.1162 26.4933 15.2837 26.1558 15.3238L25.9998 15.3332H24.6665C24.3267 15.3328 23.9998 15.2027 23.7527 14.9694C23.5056 14.7361 23.3568 14.4172 23.3369 14.078C23.317 13.7387 23.4274 13.4047 23.6455 13.1441C23.8636 12.8835 24.173 12.716 24.5105 12.6758L24.6665 12.6665H25.9998ZM6.28382 4.41317L6.40916 4.52384L7.34249 5.45717C7.58162 5.69711 7.72046 6.01908 7.7308 6.35768C7.74114 6.69628 7.62221 7.02612 7.39816 7.28021C7.17412 7.5343 6.86176 7.69358 6.52453 7.7257C6.1873 7.75783 5.85048 7.66039 5.58249 7.45317L5.45716 7.3425L4.52383 6.40917C4.29461 6.17941 4.15707 5.87394 4.13698 5.55001C4.11689 5.22609 4.21564 4.90596 4.41471 4.64964C4.61378 
        4.39332 4.8995 4.21841 5.21832 4.15769C5.53713 4.09698 5.86715 4.15463 6.14649 4.31984L6.28382 4.41317ZM23.4758 4.52384C23.7054 4.75343 23.8433 5.05889 23.8637 5.38293C23.8841 5.70696 23.7855 6.0273 23.5865 6.28384L23.4758 6.40917L22.5425 7.3425C22.3025 7.58163 21.9806 7.72047 21.642 7.73081C21.3034 7.74115 20.9735 7.62222 20.7195 7.39817C20.4654 7.17413 20.3061 6.86177 20.274 6.52454C20.2418 6.18731 20.3393 5.85049 20.5465 5.5825L20.6572 5.45717L21.5905 4.52384C21.8405 4.27388 22.1796 4.13346 22.5332 4.13346C22.8867 4.13346 23.2258 4.27388 23.4758 4.52384ZM13.9998 0.666504C14.3264 0.666547 14.6416 0.786446 14.8857 1.00346C15.1297 1.22047 15.2856 1.5195 15.3238 1.84384L15.3332 1.99984V3.33317C15.3328 3.67301 15.2027 3.99988 14.9694 4.247C14.7361 4.49411 14.4172 4.64282 14.078 4.66273C13.7387 4.68265 13.4047 4.57227 13.144 4.35415C12.8834 4.13603 12.716 3.82663 12.6758 3.48917L12.6665 3.33317V1.99984C12.6665 1.64622 12.807 1.30708 13.057 1.05703C13.3071 0.80698 13.6462 0.666504 13.9998 0.666504ZM13.9998 7.33317C15.3061 7.3331 16.5835 7.71676 17.6736 8.4365C18.7636 9.15625 19.6182 10.1804 20.1312 11.3817C20.6442 12.5829 20.793 13.9085 20.5591 15.1936C20.3252 16.4787 19.719 17.6668 18.8156 18.6103C17.9122 19.5538 16.7516 20.2112 15.4779 20.5007C14.2041 20.7902 12.8734 20.6992 11.6509 20.2389C10.4285 19.7786 9.36821 18.9693 8.60177 17.9115C7.83533 16.8538 7.39652 15.5942 7.33983 14.2892L7.33316 
        13.9998L7.33983 13.7105C7.41438 11.9942 8.14863 10.3729 9.38946 9.18472C10.6303 7.99656 12.2819 7.33327 13.9998 7.33317Z" fill="#F7CC74"/>
        </svg>`
        let dark = document.getElementById("svg");
        let darkpath = document.getElementById("pat");
        dark.style.marginTop = "3px";
        dark.style.cursor = "pointer";
        let light = document.getElementById("svgL");
        let lightpath = document.getElementById("patL");
        light.style.marginTop = "3px";
        light.style.cursor = "pointer";
        dark.onclick = () => {
            darkpath.style.fill = "#2E1A47";
            lightpath.style.fill = "#8C898E";
        }
        light.onclick = () => {
            darkpath.style.fill = "#8C898E";
            lightpath.style.fill = "#F7CC74";
        }
    }
}

function accesibility_change() {
    if (dijete.hasChildNodes() && roditelj.hasChildNodes()) {
        dijete.innerHTML = "";
        dijete.innerHTML = ` 
        <div class="outside">
            <div class="inside">
            </div>
        </div>
        <p class="pID" id="disleksija">Dyslexia</p>
        <div class="outside">
            <div class="inside">
            </div>
        </div>
        <p class="pID" id="kontrast">Contrast</p>
        <div class="outside">
            <div class="inside">
            </div>
        </div>
        <p class="pID" id="linkovi">Mark Links</p>
        `
        const gumb = document.querySelectorAll(".inside");
        gumb.forEach((el) => {
            el.onclick = () => {
                if (el.classList.contains("inside_click")) {
                    el.classList.remove("inside_click");
                    el.classList.add("reverse_inside_click");
                } else {
                    el.classList.remove("reverse_inside_click");
                    el.classList.add("inside_click");
                }
            }
        });
        const disleksija = document.getElementById("disleksija");
        const kontrast = document.getElementById("kontrast");
        const linkovi = document.getElementById("linkovi");
        if (document.documentElement.lang == "hr") {
            disleksija.textContent = "Disleksija"
            kontrast.textContent = "Kontrast"
            linkovi.textContent = "Označi linkove"
        } else if (document.documentElement.lang == "de") {

        }
    }
}

language.addEventListener("click", language_change);
theme.addEventListener("click", theme_change);
accesibility.addEventListener("click", accesibility_change);

button.addEventListener("click", () => {
    svg.classList.toggle("change");
    if (dropdown.classList.contains("show")) {
        dropdown.classList.remove("show");
        dropdown.classList.add("unshow");
        dropdownP.forEach((el) => {
            setTimeout(() => {
                el.style.display = "none";
            }, 1000);
        });
        but.forEach((element) => {
            setTimeout(() => {
                element.style.display = "none";
            }, 1000);
        });
        seperator.forEach((node) => {
            setTimeout(() => {
                node.style.display = "none";
            }, 1000);
        });
        setTimeout(() => {
            dijete.style.display = "none";
        }, 1000);

    } else {
        dropdown.classList.remove("unshow");
        dropdown.classList.add("show");
        dropdownP.forEach((el) => {
            el.style.display = "inline-block";
        });
        but.forEach((element) => {
            element.style.display = "inline-block";
        });
        seperator.forEach((node) => {
            node.style.display = "inline-block";
        });
        dijete.style.display = "flex";

    }
})

let animateFlag = true;
let animateFlagMiddle = true;
let animateFlagBottom = true;
let zastavica = true;
let element = document.getElementById("top").offsetTop;
let elementMiddle = document.getElementById("middle").offsetTop;
let elementBottom = document.getElementById("bottom").offsetTop;
let ficr = document.querySelector(".FeaturesSection")
let FeatureSection = document.getElementById("Pocetak");
let dva = document.getElementById("dva");
let tri = document.getElementById("tri");
let elementFeature = FeatureSection.offsetTop;
let elementScroll = document.getElementById("top");
let elementScrollMiddle = document.getElementById("middle");
let elementScrollBottom = document.getElementById("bottom");
let screenHeight = window.innerHeight;
let activationOffset = 0.5;
let activationPoint = element - (screenHeight * activationOffset);
let activationPointMiddle = elementMiddle - (screenHeight * activationOffset);
let activationPointBottom = elementBottom - (screenHeight * activationOffset);
let activationPointFeature = elementFeature - (screenHeight * activationOffset);

window.addEventListener("scroll", function () {
    if (this.pageYOffset > activationPoint) {
        if (animateFlag) {
            elementScroll.classList.add("MoveRight");
            animateFlag = false;
        }
    }
    if (this.pageYOffset > activationPointMiddle) {
        if (animateFlagMiddle) {
            elementScrollMiddle.classList.add("MoveLeftMiddle");
            animateFlagMiddle = false;
        }
    }
    if (this.pageYOffset > activationPointBottom) {
        if (animateFlagBottom) {
            elementScrollBottom.classList.add("MoveRight");
            animateFlagBottom = false;
        }
    }
    if (this.pageYOffset > activationPointFeature) {
        if (zastavica) {
            FeatureSection.classList.add("goUp");
            ficr.classList.add("paddingChange");
            dva.classList.add("goUp");
            tri.classList.add("goUp");
            zastavica = false;
        }
    }
})




