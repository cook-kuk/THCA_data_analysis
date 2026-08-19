import React from "react";

interface Props {
  title: string;
  body?: string;
  src?: string;
  alt?: string;
  stat?: string;
  note?: string;
  onOpen?: (src: string, alt: string) => void;
}

export const FigureCard: React.FC<Props> = ({ title, body, src, alt, stat, note, onOpen }) => (
  <div className="fig-card">
    <div className={`fig-thumb${src ? "" : " placeholder"}`}>
      {src ? (
        <img
          src={src}
          alt={alt || title}
          onClick={() => onOpen && onOpen(src, alt || title)}
        />
      ) : (
        <span>{title} placeholder — replace with final PNG</span>
      )}
    </div>
    <div className="fig-title">{title}</div>
    {body && <p className="fig-body">{body}</p>}
    {stat && <span className="fig-stat">{stat}</span>}
    {note && <div className="fig-note">{note}</div>}
  </div>
);
