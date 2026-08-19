/* THCA Hub — minimal lightbox (figure click → modal zoom) */
(function(){
  'use strict';
  function init(){
    var overlay = document.createElement('div');
    overlay.className = 'lightbox-overlay';
    overlay.innerHTML = '<span class="lightbox-close" aria-label="close">×</span><img alt=""/><div class="lightbox-caption"></div>';
    document.body.appendChild(overlay);
    var modalImg = overlay.querySelector('img');
    var modalCap = overlay.querySelector('.lightbox-caption');

    function close(){ overlay.classList.remove('show'); }
    overlay.addEventListener('click', close);
    overlay.querySelector('.lightbox-close').addEventListener('click', function(e){e.stopPropagation(); close();});
    document.addEventListener('keydown', function(e){ if(e.key === 'Escape') close(); });

    document.querySelectorAll('figure img').forEach(function(img){
      img.addEventListener('click', function(){
        modalImg.src = this.src;
        modalImg.alt = this.alt || '';
        var cap = this.closest('figure').querySelector('figcaption');
        modalCap.textContent = cap ? cap.textContent.trim().slice(0, 240) : '';
        modalCap.style.display = cap ? 'block' : 'none';
        overlay.classList.add('show');
      });
    });
  }
  if(document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
