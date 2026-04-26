/* Section 1 — Three.js 3D PCA scatter with pre/post ComBat flip animation. */
(function (global) {
  const CLASS_A = '#F5A623';   // amber
  const CLASS_B = '#6AB04C';   // green

  let scene, camera, renderer, controls;
  let pointCloudSolid = null;   // TCGA (filled discs)
  let pointCloudRing = null;    // GEO  (rings)
  let axisLine = null;
  let container = null;
  let tweenState = null;
  let rafId = null;

  function makeCircleTexture(filled) {
    const size = 64;
    const c = document.createElement('canvas');
    c.width = c.height = size;
    const ctx = c.getContext('2d');
    ctx.clearRect(0, 0, size, size);
    ctx.beginPath();
    ctx.arc(size / 2, size / 2, size / 2 - 6, 0, Math.PI * 2);
    if (filled) {
      ctx.fillStyle = '#ffffff';
      ctx.fill();
    } else {
      ctx.lineWidth = 6;
      ctx.strokeStyle = '#ffffff';
      ctx.stroke();
    }
    const tex = new THREE.CanvasTexture(c);
    tex.minFilter = THREE.LinearFilter;
    return tex;
  }

  function initScene() {
    container = document.getElementById('three-host');
    if (!container) return;
    const w = container.clientWidth;
    const h = container.clientHeight || 560;

    scene = new THREE.Scene();
    scene.background = null;

    camera = new THREE.PerspectiveCamera(55, w / h, 0.01, 100);
    camera.position.set(2.4, 1.8, 2.8);

    renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize(w, h);
    renderer.setClearColor(0x000000, 0);
    container.appendChild(renderer.domElement);

    controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.08;
    controls.rotateSpeed = 0.55;

    // ambient + grid
    scene.add(new THREE.AmbientLight(0xffffff, 0.9));

    const grid = new THREE.GridHelper(3, 6, 0x2a2a30, 0x1a1a20);
    grid.position.y = -1.05;
    scene.add(grid);

    // axis line placeholder
    const axisMat = new THREE.LineDashedMaterial({ color: 0xF5A623, dashSize: 0.08, gapSize: 0.06, transparent: true, opacity: 0.8 });
    const axisGeom = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(-1.2, 0, 0), new THREE.Vector3(1.2, 0, 0)]);
    axisLine = new THREE.Line(axisGeom, axisMat);
    axisLine.computeLineDistances();
    scene.add(axisLine);

    // resize
    const onResize = () => {
      const ww = container.clientWidth;
      const hh = container.clientHeight || 560;
      renderer.setSize(ww, hh);
      camera.aspect = ww / hh;
      camera.updateProjectionMatrix();
    };
    window.addEventListener('resize', onResize);

    // events from controls overlay
    container.addEventListener('three-reset', () => {
      controls.reset();
      camera.position.set(2.4, 1.8, 2.8);
    });
    container.addEventListener('three-tour', () => {
      const app = window.__v9;
      if (!app) return;
      app.toggleFlip();
      setTimeout(() => app.toggleFlip(), 2500);
    });

    // start loop
    (function loop() {
      rafId = requestAnimationFrame(loop);
      controls.update();
      // slow auto rotate when idle
      if (tweenState && tweenState.active) {
        const t = (performance.now() - tweenState.start) / tweenState.duration;
        const clamp = Math.max(0, Math.min(1, t));
        const eased = clamp < 0.5 ? 2 * clamp * clamp : 1 - Math.pow(-2 * clamp + 2, 2) / 2;
        const from = tweenState.from, to = tweenState.to, arr = tweenState.posAttr.array;
        for (let i = 0; i < arr.length; i++) {
          arr[i] = from[i] + (to[i] - from[i]) * eased;
        }
        tweenState.posAttr.needsUpdate = true;
        // axis rotation tween
        if (tweenState.axisFrom !== undefined) {
          const ang = tweenState.axisFrom + (tweenState.axisTo - tweenState.axisFrom) * eased;
          axisLine.rotation.z = ang;
        }
        if (clamp >= 1) tweenState.active = false;
      }
      renderer.render(scene, camera);
    })();
  }

  function buildPointClouds(cancer, pcaData) {
    if (!scene) return;
    if (!pcaData || !pcaData.pre) return;
    // remove old
    if (pointCloudSolid) { scene.remove(pointCloudSolid); pointCloudSolid.geometry.dispose(); pointCloudSolid.material.dispose(); }
    if (pointCloudRing) { scene.remove(pointCloudRing); pointCloudRing.geometry.dispose(); pointCloudRing.material.dispose(); }

    const coords = pcaData.pre || [];
    const classA = window.__v9?.cohortInfo?.class_a || 'A';

    const solidPos = [], solidCol = [];
    const ringPos = [], ringCol = [];

    for (const p of coords) {
      const [x, y, z, cohort, label] = p;
      const isTCGA = String(cohort).startsWith('TCGA');
      const col = new THREE.Color(label === classA ? CLASS_A : CLASS_B);
      if (isTCGA) {
        solidPos.push(x, y, z);
        solidCol.push(col.r, col.g, col.b);
      } else {
        ringPos.push(x, y, z);
        ringCol.push(col.r, col.g, col.b);
      }
    }

    // solid points
    if (solidPos.length) {
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.Float32BufferAttribute(solidPos, 3));
      g.setAttribute('color', new THREE.Float32BufferAttribute(solidCol, 3));
      const tex = makeCircleTexture(true);
      const m = new THREE.PointsMaterial({
        size: 0.09,
        map: tex, alphaTest: 0.3, transparent: true,
        vertexColors: true,
        sizeAttenuation: true,
      });
      pointCloudSolid = new THREE.Points(g, m);
      scene.add(pointCloudSolid);
    }
    if (ringPos.length) {
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.Float32BufferAttribute(ringPos, 3));
      g.setAttribute('color', new THREE.Float32BufferAttribute(ringCol, 3));
      const tex = makeCircleTexture(false);
      const m = new THREE.PointsMaterial({
        size: 0.11,
        map: tex, alphaTest: 0.3, transparent: true,
        vertexColors: true,
        sizeAttenuation: true,
      });
      pointCloudRing = new THREE.Points(g, m);
      scene.add(pointCloudRing);
    }
  }

  function tweenToMode(combatOn, pcaData) {
    if (!scene || !axisLine) return;  // scene not initialized yet
    if (!pointCloudSolid && !pointCloudRing) return;
    const coords = combatOn ? (pcaData.post || []) : (pcaData.pre || []);
    const classA = window.__v9?.cohortInfo?.class_a || 'A';
    const targetSolid = [], targetRing = [];
    for (const p of coords) {
      const [x, y, z, cohort] = p;
      const isTCGA = String(cohort).startsWith('TCGA');
      if (isTCGA) targetSolid.push(x, y, z);
      else targetRing.push(x, y, z);
    }

    // tween both clouds concurrently by stashing their separate arrays
    const start = performance.now();
    const duration = 2000;

    // solid
    if (pointCloudSolid) {
      const posAttr = pointCloudSolid.geometry.attributes.position;
      const fromArr = Array.from(posAttr.array);
      // pad/trim if count differs (shouldn't, but safety)
      const toArr = new Array(fromArr.length);
      for (let i = 0; i < fromArr.length; i++) toArr[i] = targetSolid[i] != null ? targetSolid[i] : fromArr[i];
      animateArray(posAttr, fromArr, toArr, start, duration);
    }
    if (pointCloudRing) {
      const posAttr = pointCloudRing.geometry.attributes.position;
      const fromArr = Array.from(posAttr.array);
      const toArr = new Array(fromArr.length);
      for (let i = 0; i < fromArr.length; i++) toArr[i] = targetRing[i] != null ? targetRing[i] : fromArr[i];
      animateArray(posAttr, fromArr, toArr, start, duration);
    }

    // axis rotation — pi/2 flip for visual "inversion"
    const axFrom = axisLine.rotation.z;
    const axTo = combatOn ? Math.PI / 2 : 0;
    (function ax() {
      const t = Math.min(1, (performance.now() - start) / duration);
      const eased = t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;
      axisLine.rotation.z = axFrom + (axTo - axFrom) * eased;
      if (t < 1) requestAnimationFrame(ax);
    })();
  }

  function animateArray(posAttr, fromArr, toArr, start, duration) {
    (function step() {
      const t = Math.min(1, (performance.now() - start) / duration);
      const eased = t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;
      const a = posAttr.array;
      for (let i = 0; i < a.length; i++) a[i] = fromArr[i] + (toArr[i] - fromArr[i]) * eased;
      posAttr.needsUpdate = true;
      if (t < 1) requestAnimationFrame(step);
    })();
  }

  const Section1 = {
    init(app) {
      // Wait one tick so Alpine renders the container
      setTimeout(() => {
        initScene();
        const pca = app.selectedPCA;
        buildPointClouds(app.selectedCancer, pca);
        tweenToMode(app.combatOn, pca);
      }, 60);
    },
    onSelection(app) {
      const pca = app.selectedPCA;
      buildPointClouds(app.selectedCancer, pca);
      tweenToMode(app.combatOn, pca);
    },
    onFlip(app) {
      tweenToMode(app.combatOn, app.selectedPCA);
    },
  };

  global.Section1 = Section1;
})(window);
