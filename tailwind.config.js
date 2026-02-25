/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  darkMode: 'class',
  theme: {
    extend: {
      colors: {
        primary: {
          50: '#f0f7ff',
          100: '#e0efff',
          200: '#b8dbff',
          300: '#7abfff',
          400: '#3a9fff',
          500: '#0a7fff',
          600: '#0062d6',
          700: '#004dad',
          800: '#004190',
          900: '#003876',
        },
        fossil: {
          light: '#f5e6d3',
          DEFAULT: '#c4956a',
          dark: '#8b6040',
        },
      },
    },
  },
  plugins: [],
}

